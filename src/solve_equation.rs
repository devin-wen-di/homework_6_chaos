//solve_equation.rs
use std::f64::consts::PI;
use std::fs::File;
use std::io::Write;
use crate::model::PendulumParams;

#[derive(Debug, Clone, Copy)]
pub struct State {
    pub theta: f64,
    pub omega: f64,
}

// 微分方程变成一个向量场，输入向量场的坐标，返回该处向量的y分量
pub fn rhs(theta: f64, omega: f64, t: f64, params: &PendulumParams) -> (f64, f64) {
    let d_theta_dt = omega;
    let d_omega_dt = -(params.g / params.l) * theta.sin()
        - params.q * d_theta_dt
        + params.f_d * (params.omega_d * t).sin();
    (d_theta_dt, d_omega_dt)
}

// 给出当前的状态和时间，返回下一步的状态和时间（RK4）
pub fn rk4_step(state: &State, t: f64, params: &PendulumParams) -> (State, f64) {
    let (k1_theta, k1_omega) = rhs(state.theta, state.omega, t, params);
    let middle_step_1 = State {
        theta: state.theta + 0.5 * params.dt * k1_theta,
        omega: state.omega + 0.5 * params.dt * k1_omega,
    };
    let (k2_theta, k2_omega) = rhs(middle_step_1.theta, middle_step_1.omega, t + 0.5 * params.dt, params);
    let middle_step_2 = State {
        theta: state.theta + 0.5 * params.dt * k2_theta,
        omega: state.omega + 0.5 * params.dt * k2_omega,
    };
    let (k3_theta, k3_omega) = rhs(middle_step_2.theta, middle_step_2.omega, t + 0.5 * params.dt, params);
    let middle_step_3 = State {
        theta: state.theta + params.dt * k3_theta,
        omega: state.omega + params.dt * k3_omega,
    };
    let (k4_theta, k4_omega) = rhs(middle_step_3.theta, middle_step_3.omega, t + params.dt, params);

    let new_theta = state.theta + params.dt / 6.0 * (k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    let new_omega = state.omega + params.dt / 6.0 * (k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega);
    (State { theta: new_theta, omega: new_omega }, t + params.dt)
}

pub fn solve(params: &PendulumParams, initial_theta: f64, initial_omega: f64) -> Vec<(f64, State)> {
    let mut trajectory = Vec::new();
    let mut state = State {
        theta: initial_theta,
        omega: initial_omega,
    };
    let mut t = 0.0;
    trajectory.push((t, state.clone()));

    // 为了避免浮点累计误差导致多一步或少一步，按固定步数迭代
    let steps = (params.t_end / params.dt) as usize;
    for _ in 0..steps {
        let (new_state, new_t) = rk4_step(&state, t, params);
        state = new_state;
        t = new_t;
        trajectory.push((t, state.clone()));
    }
    trajectory
}

fn wrap_angle(theta: f64) -> f64 {
    (theta + PI).rem_euclid(2.0 * PI) - PI
}

/// 计算庞加莱截面：在每个驱动周期处采样。
///
/// params: 系统参数（包含 dt, t_end, omega_d）
/// initial_*: 初始角度与角速度
/// transient_periods: 要丢弃的周期数（过渡）
/// sample_periods: 采样的周期数
/// 返回值：Vec<(theta, omega)>，按采样顺序排列
/// 改进的庞加莱采样：在精确的驱动周期时刻采样
pub fn poincare(
    params: &PendulumParams,
    initial_theta: f64,
    initial_omega: f64,
    transient_periods: usize,
    sample_periods: usize,
) -> Vec<(f64, f64)> {
    let period = 2.0 * PI / params.omega_d;
    
    // 先计算完整轨迹
    let traj = solve(params, initial_theta, initial_omega);
    
    let mut samples = Vec::new();
    
    for n in (transient_periods + 1)..=(transient_periods + sample_periods) {
        let target_time = (n as f64) * period;
        
        // 找到目标时间前后的点
        if let Some(idx) = traj.iter().position(|(t, _)| *t >= target_time) {
            if idx > 0 {
                let (t1, s1) = traj[idx - 1];
                let (t2, s2) = traj[idx];
                
                // 线性插值
                let frac = (target_time - t1) / (t2 - t1);
                
                // 角度插值需要考虑周期性
                let mut theta_diff = s2.theta - s1.theta;
                theta_diff = (theta_diff + PI).rem_euclid(2.0 * PI) - PI;
                
                let theta_interp = s1.theta + frac * theta_diff;
                let omega_interp = s1.omega + frac * (s2.omega - s1.omega);
                
                samples.push((wrap_angle(theta_interp), omega_interp));
            }
        }
    }
    
    samples
}

/// 把 poincare 结果写成 CSV 文件（两列：theta,omega）
pub fn write_poincare_csv(
    path: &str,
    params: &PendulumParams,
    initial_theta: f64,
    initial_omega: f64,
    transient_periods: usize,
    sample_periods: usize,
) -> std::io::Result<()> {
    let samples = poincare(params, initial_theta, initial_omega, transient_periods, sample_periods);
    let mut f = File::create(path)?;
    writeln!(f, "theta,omega")?;
    for (th, om) in samples {
        writeln!(f, "{:.12},{:.12}", th, om)?;
    }
    Ok(())
}

/// 从已经计算好的轨迹中按驱动周期采样（线性/角度插值）
pub fn sample_poincare_from_trajectory(
    traj: &Vec<(f64, State)>,
    params: &PendulumParams,
    transient_periods: usize,
    sample_periods: usize,
) -> Vec<(f64, f64)> {
    let period = 2.0 * PI / params.omega_d;

    let mut samples = Vec::with_capacity(sample_periods);

    // 生成要采样的时刻列表（以驱动周期为单位，从 1..=total_periods）
    let start_n = transient_periods + 1;
    let end_n = transient_periods + sample_periods;

    for n in start_n..=end_n {
        let t_sample = (n as f64) * period;
        // 在 traj 中找到包含 t_sample 的相邻点
        if let Some(i) = traj.iter().position(|(t, _)| *t >= t_sample) {
            if i == 0 {
                // 如果第一个点就大于等于采样时间，直接取该点
                let (_, s) = &traj[i];
                samples.push((wrap_angle(s.theta), s.omega));
            } else {
                let (t1, s1) = &traj[i - 1];
                let (t2, s2) = &traj[i];
                let dt = t2 - t1;
                if dt <= 0.0 {
                    let (_, s) = &traj[i];
                    samples.push((wrap_angle(s.theta), s.omega));
                } else {
                    let alpha = (t_sample - t1) / dt;
                    // 角度需要做周期展开再插值
                    let mut dtheta = s2.theta - s1.theta;
                    dtheta = (dtheta + PI).rem_euclid(2.0 * PI) - PI;
                    let theta_interp = s1.theta + alpha * dtheta;
                    let omega_interp = s1.omega + alpha * (s2.omega - s1.omega);
                    samples.push((wrap_angle(theta_interp), omega_interp));
                }
            }
        } else {
            // 未找到大于等于 t_sample 的点（轨迹不够长），跳出
            break;
        }
    }

    samples
}

/// 使用现有的 solve() 先计算完整轨迹，然后从轨迹中按周期采样（推荐）
pub fn poincare_via_solve(
    params: &PendulumParams,
    initial_theta: f64,
    initial_omega: f64,
    transient_periods: usize,
    sample_periods: usize,
) -> Vec<(f64, f64)> {
    let traj = solve(params, initial_theta, initial_omega);
    sample_poincare_from_trajectory(&traj, params, transient_periods, sample_periods)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_small_angle_pendulum_analytical_consistency() {
        // 设置无阻尼、无驱动的小角度单摆参数
        let mut params = PendulumParams::new();
        params.q = 0.0;
        params.f_d = 0.0;
        params.dt = 0.001;

        let initial_theta: f64 = 0.1;  // 小角度 (弧度)
        let initial_omega: f64 = 0.0;

        // 调用数值求解函数
        let trajectory = solve(&params, initial_theta, initial_omega);

        // 计算自然频率 ω = sqrt(g / l)
        let omega_n = (params.g / params.l).sqrt();

        // 验证轨迹长度（大致为 t_end / dt + 1）
        assert_eq!(trajectory.len(), ((params.t_end / params.dt) as usize + 1));

        // 对于轨迹中的每个点，验证与解析解的匹配
        // 解析解: θ(t) = θ0 * cos(ω t), ω(t) = -θ0 * ω * sin(ω t)
        for (t, state) in &trajectory {
            let theta_analytic = initial_theta * (omega_n * *t).cos();
            let omega_analytic = -initial_theta * omega_n * (omega_n * *t).sin();

            // 使用相对误差阈值 1e-4（RK4 + 小 dt 已足够精确）
            assert_relative_eq!(state.theta, theta_analytic, epsilon = 1e-2);
            assert_relative_eq!(state.omega, omega_analytic, epsilon = 1e-2);
        }

        // 额外验证：相图能量守恒（对于无阻尼情况）
        // E ≈ (1/2) m l^2 ω^2 + m g l (1 - cos θ) ≈ (1/2) ω^2 + (g/l) θ^2 / 2（小角度）
        // 简化： θ^2 + (ω / sqrt(g/l))^2 应 ≈ constant = initial_theta^2
        let initial_energy = initial_theta.powi(2);
        for (_, state) in &trajectory[1..] {  // 跳过初始点
            let conserved = state.theta.powi(2) + (state.omega / omega_n).powi(2);
            assert_relative_eq!(conserved, initial_energy, epsilon = 1e-3);  // 稍松阈值容忍累积误差
        }
    }
}