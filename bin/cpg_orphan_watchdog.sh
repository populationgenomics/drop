#!/usr/bin/env bash
# CPG: kill whatever a task's main script leaves running once it has exited (docs/CPG_METAMIST_SEQR_DESIGN.md, T5.7).
#
# Started in the background by the beforeScript in conf/cpg.config, detached from the task's stdin/stdout/stderr.
# Without it, a forked worker that outlives its parent (e.g. BiocParallel workers after one is OOM-killed) keeps the
# wrapper's `nxf_launch | tee` pipes open, so `wait $pid` in .command.run blocks and the task and its VM stay up.
#
#   cpg_orphan_watchdog.sh <task work dir> [log file]
#
# Finds the main process (the earliest one whose command line contains <task work dir>/.command.sh and whose parent's
# doesn't) and records its descendants with their /proc start times while it runs. Once it has exited, a survivor is
# never killed while it is writing output: local disk writes (write_bytes rising), writes to a file it holds open
# for writing on a cloud mount, or writes over a socket (wchar rising). Only survivors idle for CPG_WATCHDOG_GRACE
# seconds get TERM, then KILL. If a process's /proc/<pid>/io can't be read, it is treated as writing.
# Uses only bash and /proc: no ps/pkill in these images.

task_dir="${1:?usage: cpg_orphan_watchdog.sh <task work dir> [log file]}"
log="${2:-/dev/null}"
interval="${CPG_WATCHDOG_INTERVAL:-10}"
grace="${CPG_WATCHDOG_GRACE:-300}"
read -r -a mounts <<< "${CPG_WATCHDOG_MOUNTS:-/mnt/disks/ /fusion/}"
# inputs are staged after the beforeScript starts us, and copying a large BAM (T5.5) can take a while
startup_timeout="${CPG_WATCHDOG_STARTUP_TIMEOUT:-3600}"
main_script="${task_dir%/}/.command.sh"
self=$$

# no /proc (e.g. macOS): nothing to watch
[[ -r /proc/self/stat ]] || exit 0
command -v readlink > /dev/null 2>&1 && HAVE_READLINK=1 || HAVE_READLINK=0

say() {
    printf '%s cpg_orphan_watchdog: %s\n' "$(date -u +%FT%TZ)" "$*" >> "$log" 2>/dev/null || true
}

# sets STAT_PPID and STAT_START for a pid; fails if it has gone
read_stat() {
    local line fields
    read -r line < "/proc/$1/stat" 2>/dev/null || return 1
    # the command name is in parentheses and may contain spaces: take the fields after its last ')'
    read -r -a fields <<< "${line##*) }"
    STAT_PPID=${fields[1]}
    STAT_START=${fields[19]}
}

cmdline_of() {
    local part cmd=""
    [[ -r /proc/$1/cmdline ]] || return 1
    while IFS= read -r -d '' part; do cmd+="$part "; done < "/proc/$1/cmdline" 2>/dev/null
    printf '%s' "$cmd"
}

is_alive() {  # pid, start time
    read_stat "$1" && [[ $STAT_START == "$2" ]]
}

main=""
main_start=""
find_main() {
    local dir pid cmd parent_cmd best_start=""
    for dir in /proc/[0-9]*; do
        pid=${dir#/proc/}
        [[ $pid == "$self" ]] && continue
        cmd=$(cmdline_of "$pid") || continue
        [[ $cmd == *"$main_script"* ]] || continue
        read_stat "$pid" || continue
        local start=$STAT_START ppid=$STAT_PPID
        parent_cmd=$(cmdline_of "$ppid" 2>/dev/null || true)
        # a forked worker shares its parent's command line; the main process is the one whose parent doesn't
        [[ $parent_cmd == *"$main_script"* ]] && continue
        if [[ -z $best_start || $start -lt $best_start ]]; then
            main=$pid
            main_start=$start
            best_start=$start
        fi
    done
}

declare -A seen=()  # pid -> start time, every descendant of the main process seen so far

# add every current descendant of the given pids to `seen`
record_descendants() {
    local -A children=() starts=()
    local dir pid queue=() current child
    for dir in /proc/[0-9]*; do
        pid=${dir#/proc/}
        read_stat "$pid" || continue
        children[$STAT_PPID]+=" $pid"
        starts[$pid]=$STAT_START
    done
    queue=("$@")
    while ((${#queue[@]})); do
        current=${queue[0]}
        queue=("${queue[@]:1}")
        for child in ${children[$current]:-}; do
            [[ $child == "$self" ]] && continue
            [[ -n ${seen[$child]:-} ]] && continue
            seen[$child]=${starts[$child]}
            queue+=("$child")
        done
    done
}

survivors() {
    local pid
    for pid in "${!seen[@]}"; do
        is_alive "$pid" "${seen[$pid]}" && printf '%s ' "$pid"
    done
}

deadline=$((SECONDS + startup_timeout))
until find_main && [[ -n $main ]]; do
    if ((SECONDS > deadline)); then
        say "no process running $main_script after ${startup_timeout}s; exiting"
        exit 0
    fi
    sleep 1
done
say "watching main process $main ($(cmdline_of "$main"))"

while is_alive "$main" "$main_start"; do
    record_descendants "$main"
    sleep "$interval"
done

# workers started after the last pass, or by orphans since, are descendants of those already seen
record_descendants "${!seen[@]}"
left=$(survivors)
if [[ -z $left ]]; then
    say "main process $main exited; nothing left running"
    exit 0
fi

say "main process $main exited; still running: ${left% }; stopping any that stay idle for ${grace}s"

declare -A prev_wchar=() prev_wbytes=() idle_since=() writing_logged=()

# sets IO_WCHAR and IO_WBYTES; fails if /proc/<pid>/io can't be read
read_io() {
    local key value
    IO_WCHAR=""
    IO_WBYTES=""
    while read -r key value; do
        case $key in
            wchar:) IO_WCHAR=$value ;;
            write_bytes:) IO_WBYTES=$value ;;
        esac
    done < "/proc/$1/io" 2>/dev/null
    [[ -n $IO_WCHAR && -n $IO_WBYTES ]]
}

# sets FD_TCP and FD_MOUNT_WRITE (1/0): an established TCP connection held (possibly an upload; send() calls
# don't show up in wchar, so the connection itself is the signal), a file on a cloud mount open for writing
read_fds() {
    local fd target flags mount inode key value inodes=" "
    FD_TCP=0
    FD_MOUNT_WRITE=0
    # without readlink the fds can't be inspected: assume the process may be writing
    if ((!HAVE_READLINK)); then
        FD_TCP=1
        return 0
    fi
    for fd in /proc/"$1"/fd/*; do
        target=$(readlink "$fd" 2>/dev/null) || continue
        if [[ $target == socket:\[*\] ]]; then
            inode=${target#socket:[}
            inodes+="${inode%]} "
            continue
        fi
        for mount in "${mounts[@]}"; do
            [[ $target == "$mount"* ]] || continue
            flags=""
            while read -r key value; do
                [[ $key == flags: ]] && flags=$value
            done < "/proc/$1/fdinfo/${fd##*/}" 2>/dev/null
            # O_WRONLY = 1, O_RDWR = 2
            [[ -n $flags ]] && (( (8#$flags & 3) != 0 )) && FD_MOUNT_WRITE=1
        done
    done
    [[ $inodes == " " ]] && return 0
    # /proc/<pid>/net/tcp{,6}: field 4 is the state (01 = ESTABLISHED), field 10 the socket inode
    local table line fields
    for table in /proc/"$1"/net/tcp /proc/"$1"/net/tcp6; do
        [[ -r $table ]] || continue
        {
            read -r line  # header
            while read -r line; do
                read -r -a fields <<< "$line"
                [[ ${fields[3]:-} == 01 && $inodes == *" ${fields[9]:-x} "* ]] && { FD_TCP=1; return 0; }
            done
        } < "$table" 2>/dev/null
    done
}

# 0 if the pid wrote output since the last check (or can't be checked), with the reason in WRITING_REASON;
# records the counters for next time
is_writing() {
    local pid=$1 before_wchar=${prev_wchar[$1]:-} before_wbytes=${prev_wbytes[$1]:-}
    WRITING_REASON=""
    if ! read_io "$pid"; then
        WRITING_REASON="its /proc/$pid/io can't be read"
        return 0
    fi
    prev_wchar[$pid]=$IO_WCHAR
    prev_wbytes[$pid]=$IO_WBYTES
    # the first look only sets the baseline
    if [[ -z $before_wchar ]]; then
        WRITING_REASON="baseline"
        return 0
    fi
    if ((IO_WBYTES > before_wbytes)); then
        WRITING_REASON="writing to disk"
        return 0
    fi
    read_fds "$pid"
    if ((FD_MOUNT_WRITE && IO_WCHAR > before_wchar)); then
        WRITING_REASON="writing to a cloud mount"
        return 0
    fi
    if ((FD_TCP)); then
        WRITING_REASON="holding an established network connection (possibly uploading)"
        return 0
    fi
    return 1
}

while true; do
    record_descendants "${!seen[@]}"
    left=$(survivors)
    if [[ -z $left ]]; then
        say "nothing left running; done"
        exit 0
    fi
    to_stop=""
    for pid in $left; do
        if is_writing "$pid"; then
            unset "idle_since[$pid]"
            if [[ $WRITING_REASON != baseline && ${writing_logged[$pid]:-} != "$WRITING_REASON" ]]; then
                say "process $pid: $WRITING_REASON; leaving it running"
                writing_logged[$pid]=$WRITING_REASON
            fi
        else
            [[ -z ${idle_since[$pid]:-} ]] && idle_since[$pid]=$SECONDS
            ((SECONDS - idle_since[$pid] >= grace)) && to_stop+="$pid "
        fi
    done
    if [[ -n $to_stop ]]; then
        say "idle for ${grace}s, sending TERM to: ${to_stop% }"
        for pid in $to_stop; do kill -TERM "$pid" 2>/dev/null; done
        sleep 5
        for pid in $to_stop; do
            if is_alive "$pid" "${seen[$pid]}"; then
                say "sending KILL to: $pid"
                kill -KILL "$pid" 2>/dev/null
            fi
        done
    fi
    sleep "$interval"
done