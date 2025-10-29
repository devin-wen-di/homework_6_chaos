use std::f64::consts::PI;
use crate::model::PendulumParams;

#[derive(Debug, Clone, Copy)]
pub struct State{
    pub theta: f64,
    pub omega: f64,
}

//微分方程变成一个向量场，输入向量场的坐标，返回该处向量的y分量(t分量保持为1，在dt中已经体现了)
pub fn rhs(theta:f64, omega:f64, t:f64, params: $PendulumParams) -> (f64,f64) {
    let d_theta_dt = omega;
    let d_omega_dt = -(params.g/params.l)*theta.sin()
                    - params.q*d_theta_dt
                    +params.f_d*(params.omega_d*t).sin();
    (d_theta_dt, d_omega_dt)
}

//给出当前的状态和时间，返回下一步的状态和时间
pub fn rk4_step(state: &State, t: f64, params: &PendulumParams) -> (State, f64) {
    let (k1_theta,k1_omega) = rhs(state.theta, state.omega, t, params);
    let middle_stpe_1 = State {
        theta: state.theta + 0.5 * params.dt *k1_theta,
        omega: state.omega + 0.5 * params.dt *k1_omega,
    };
    let (k2_theta, k2_omega) = rhs(middle_stpe_1.theta, middle_stpe_1.omega, t + 0.5 * params.dt, params);
    let middle_stpe_2 = State {
        theta: state.theta + 0.5 * params.dt *k2_theta,
        omega: state.omega + 0.5 * params.dt *k2_omega,
    };
    let (k3_theta, k3_omega) = rhs(middle_stpe_2.theta, middle_stpe_2.omega, t + 0.5 * params.dt, params);
    let middle_stpe_3 = State {
        theta: state.theta + params.dt *k3_theta,
        omega: state.omega + params.dt *k3_omega,
    };
    let (k4_theta, k4_omega) = rhs(middle_stpe_3.theta, middle_stpe_3.omega, t + params.dt, params);

    let new_theta = state.theta + params.dt / 6.0 *(k1_theta + 2.0 * k2_theta + 2.0 * k3_theta + k4_theta);
    let new_omega = state.omega + params.dt / 6.0 *(k1_omega + 2.0 * k2_omega + 2.0 * k3_omega + k4_omega);
    (State {theta: new_theta, omega: new_omega}, t+params.dt)
}

pub fn solve(params: &PendulumParams, initial_theta: f64, initial_omega: f64) -> Vec<(f64, State)> {
    let mut trajectory = Vec::new();
    let mut state = State {
        theta: initial_theta,
        omega: initial_omega,
    };
    let mut t = 0.0;
    trajectory.push((t, state.clone()));

    while t < params.t_end {
        let (new_state, new_t) = rk4_step(state, t, params)
        state = new_state;
        t = new_t;
        trajectory.push((t, state.clone()));
    }
    trajectory
}