//model.rs
#[derive(Debug, Clone)]

pub struct PendulumParams {
    //运动方程参数
    pub g: f64,
    pub l: f64,
    pub q: f64,
    pub f_d: f64,
    pub omega_d: f64,
    
    //积分参数
    pub dt: f64,
    pub t_end: f64,

    //遍历参数
    pub theta_start: f64,
    pub theta_end: f64,
    pub d_theta: f64,
    pub omega_start: f64,
    pub omega_end: f64,
    pub d_omega: f64,
}

impl PendulumParams {
    pub fn new() -> Self {
        Self {
            g: 9.8,
            l: 1.0,
            q: 0.1,
            f_d: 1.0,
            omega_d: 1.0,

            dt: 0.01,
            t_end: 10.0,
            
            theta_start: -4.0,
            theta_end: 4.0,
            d_theta: 0.01,
            omega_start: -4.0,
            omega_end: 4.0,
            d_omega: 0.01,
        }
    }
}