use std::env;
use std::mem::size_of;
use std::os::raw::{c_int, c_long, c_uint, c_ulonglong};
use std::time::Instant;

pub(crate) struct MemoryProfiler {
    enabled: bool,
    start: Instant,
}

impl MemoryProfiler {
    pub(crate) fn from_env() -> Self {
        MemoryProfiler {
            enabled: memory_profile_enabled(),
            start: Instant::now(),
        }
    }

    pub(crate) fn checkpoint(&self, label: &str, details: impl AsRef<str>) {
        if !self.enabled {
            return;
        }

        let current_rss = current_rss_bytes()
            .map(|bytes| format!("rss={}", format_bytes(bytes)))
            .unwrap_or_else(|| "rss=unknown".to_string());
        let peak_rss = peak_rss_bytes()
            .map(|bytes| format!("peak_rss={}", format_bytes(bytes)))
            .unwrap_or_else(|| "peak_rss=unknown".to_string());
        let elapsed = self.start.elapsed().as_secs_f64();
        let details = details.as_ref();
        if details.is_empty() {
            eprintln!("[mem-profile] {elapsed:.3}s {label}: {current_rss}; {peak_rss}");
        } else {
            eprintln!("[mem-profile] {elapsed:.3}s {label}: {current_rss}; {peak_rss}; {details}");
        }
    }
}

fn memory_profile_enabled() -> bool {
    env::var("PAN_NO_REC_MEM_PROFILE")
        .map(|value| {
            let value = value.trim();
            !(value.is_empty()
                || value == "0"
                || value.eq_ignore_ascii_case("false")
                || value.eq_ignore_ascii_case("no"))
        })
        .unwrap_or(false)
}

#[repr(C)]
struct TimeVal {
    tv_sec: c_long,
    tv_usec: c_long,
}

#[repr(C)]
struct RUsage {
    ru_utime: TimeVal,
    ru_stime: TimeVal,
    ru_maxrss: c_long,
    ru_ixrss: c_long,
    ru_idrss: c_long,
    ru_isrss: c_long,
    ru_minflt: c_long,
    ru_majflt: c_long,
    ru_nswap: c_long,
    ru_inblock: c_long,
    ru_oublock: c_long,
    ru_msgsnd: c_long,
    ru_msgrcv: c_long,
    ru_nsignals: c_long,
    ru_nvcsw: c_long,
    ru_nivcsw: c_long,
}

unsafe extern "C" {
    fn getrusage(who: c_int, usage: *mut RUsage) -> c_int;
}

fn peak_rss_bytes() -> Option<usize> {
    const RUSAGE_SELF: c_int = 0;
    let mut usage = RUsage {
        ru_utime: TimeVal {
            tv_sec: 0,
            tv_usec: 0,
        },
        ru_stime: TimeVal {
            tv_sec: 0,
            tv_usec: 0,
        },
        ru_maxrss: 0,
        ru_ixrss: 0,
        ru_idrss: 0,
        ru_isrss: 0,
        ru_minflt: 0,
        ru_majflt: 0,
        ru_nswap: 0,
        ru_inblock: 0,
        ru_oublock: 0,
        ru_msgsnd: 0,
        ru_msgrcv: 0,
        ru_nsignals: 0,
        ru_nvcsw: 0,
        ru_nivcsw: 0,
    };

    let result = unsafe { getrusage(RUSAGE_SELF, &mut usage) };
    if result != 0 || usage.ru_maxrss < 0 {
        return None;
    }

    let max_rss = usage.ru_maxrss as usize;
    Some(platform_rss_units_to_bytes(max_rss))
}

#[cfg(target_os = "macos")]
#[repr(C)]
struct MachTimeValue {
    seconds: c_int,
    microseconds: c_int,
}

#[cfg(target_os = "macos")]
#[repr(C)]
struct MachTaskBasicInfo {
    virtual_size: c_ulonglong,
    resident_size: c_ulonglong,
    resident_size_max: c_ulonglong,
    user_time: MachTimeValue,
    system_time: MachTimeValue,
    policy: c_int,
    suspend_count: c_int,
}

#[cfg(target_os = "macos")]
unsafe extern "C" {
    fn mach_task_self() -> c_uint;
    fn task_info(
        target_task: c_uint,
        flavor: c_int,
        task_info_out: *mut c_int,
        task_info_outCnt: *mut c_uint,
    ) -> c_int;
}

#[cfg(target_os = "macos")]
fn current_rss_bytes() -> Option<usize> {
    const KERN_SUCCESS: c_int = 0;
    const MACH_TASK_BASIC_INFO: c_int = 20;

    let mut info = MachTaskBasicInfo {
        virtual_size: 0,
        resident_size: 0,
        resident_size_max: 0,
        user_time: MachTimeValue {
            seconds: 0,
            microseconds: 0,
        },
        system_time: MachTimeValue {
            seconds: 0,
            microseconds: 0,
        },
        policy: 0,
        suspend_count: 0,
    };
    let mut count = (size_of::<MachTaskBasicInfo>() / size_of::<c_int>()) as c_uint;

    let result = unsafe {
        task_info(
            mach_task_self(),
            MACH_TASK_BASIC_INFO,
            &mut info as *mut MachTaskBasicInfo as *mut c_int,
            &mut count,
        )
    };
    if result == KERN_SUCCESS {
        Some(info.resident_size as usize)
    } else {
        None
    }
}

#[cfg(not(target_os = "macos"))]
fn current_rss_bytes() -> Option<usize> {
    None
}

#[cfg(target_os = "macos")]
fn platform_rss_units_to_bytes(max_rss: usize) -> usize {
    max_rss
}

#[cfg(not(target_os = "macos"))]
fn platform_rss_units_to_bytes(max_rss: usize) -> usize {
    max_rss * 1024
}

pub(crate) fn format_bytes(bytes: usize) -> String {
    const KIB: f64 = 1024.0;
    const MIB: f64 = KIB * 1024.0;
    const GIB: f64 = MIB * 1024.0;

    let bytes_f = bytes as f64;
    if bytes_f >= GIB {
        format!("{:.2} GiB", bytes_f / GIB)
    } else if bytes_f >= MIB {
        format!("{:.2} MiB", bytes_f / MIB)
    } else if bytes_f >= KIB {
        format!("{:.2} KiB", bytes_f / KIB)
    } else {
        format!("{bytes} B")
    }
}
