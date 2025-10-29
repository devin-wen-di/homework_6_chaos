//main.rs
mod model;
mod solve_equation;

use crate::model::PendulumParams;
use crate::solve_equation::{poincare_via_solve, write_poincare_csv};

fn main() {
    // 设置参数与初始条件（可以根据需要调整）
    let mut params = PendulumParams::new();
    // 例如启用驱动/阻尼为默认值
    params.q = 0.5;
    params.l = 9.8;
    params.f_d = 1.2;
    params.omega_d = 2.0/3.0;
    params.t_end = 10000.0;

    let initial_theta = 1.0;
    let initial_omega = 0.0;

    // 采样设置：丢弃的过渡周期与采样周期数
    let transient_periods = 100_usize;
    let sample_periods = 2000_usize;
    let period = 2.0 * std::f64::consts::PI / params.omega_d;
    params.dt = period / 400.0;
    params.t_end = period * (transient_periods + sample_periods) as f64 + params.dt;

    // 方式 A：直接计算并写出 CSV 到 chaos/data 目录
    let out_dir = "data";
    // 确保目录存在
    if let Err(e) = std::fs::create_dir_all(out_dir) {
        eprintln!("Failed to create output directory {}: {}", out_dir, e);
        return;
    }
    let out_path = format!("{}/poincare.csv", out_dir);
    match write_poincare_csv(
        &out_path,
        &params,
        initial_theta,
        initial_omega,
        transient_periods,
        sample_periods,
    ) {
        Ok(()) => println!("Wrote Poincaré data to {}", out_path),
        Err(e) => eprintln!("Failed to write CSV: {}", e),
    }

}

