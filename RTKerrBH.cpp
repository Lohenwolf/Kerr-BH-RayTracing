#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <algorithm>

// ================= PARAMETERS ================
const double M = 1.0; // Mass
const double a = 0.999; // Spin (change this to see the BH shadow change shape) (HARD LIMIT IS 1!)
const double R_HORIZON = M + std::sqrt(M*M - a*a);

// ================= CAMERA CONFIGURATION ================
const int W = 1200;
const int H = 675;
const double FOV = 22.0; 
const double R_CAM = 1000.0;
const double THETA_CAM = 1.57079632679; // Equatorial plane

// ================= RK4 INTEGRATION SETTINGS ==================
const double MAX_DIST = 2000.0;
const int MAX_STEPS = 2000000; 

// =============== STATE ==================
struct State {
    double r, theta, phi, t;
    double pr, ptheta, pphi, pt;
};

// =============== DERIVATIVES (Boyer-Lindquist) =================
void get_derivatives(const State& s, State& d) {
    double r = s.r;
    double th = s.theta;
    double r2 = r*r;
    double a2 = a*a;

    double cost = std::cos(th);
    double sint = std::sin(th);
    double sin2 = sint*sint;
    
    // Clamp to prevent singularity at poles
    if (sin2 < 1e-9) sin2 = 1e-9;

    double Sigma = r2 + a2 * cost*cost;
    double Delta = r2 - 2.0*M*r + a2;
    if (Delta < 1e-9) Delta = 1e-9;

    double inv_grr = Delta / Sigma;
    double inv_gth = 1.0 / Sigma;

    bool zero_L = (std::abs(s.pphi) < 1e-9);
    double pphi2 = zero_L ? 0.0 : s.pphi * s.pphi;

    double inv_gtt = - ((r2 + a2)*(r2 + a2)/Delta - a2*sin2) / Sigma;
    double inv_gtp = - (a*(r2 + a2)/Delta - a) / Sigma;
    
    double term_phi = 0.0;
    double inv_gpp = 0.0;
    
    if (!zero_L) {
        inv_gpp = (1.0/sin2 - a2/Delta) / Sigma;
        term_phi = inv_gpp * s.pphi;
    }

    d.r     = inv_grr * s.pr;
    d.theta = inv_gth * s.ptheta;
    d.phi   = inv_gtp * s.pt + term_phi;
    d.t     = inv_gtt * s.pt + inv_gtp * s.pphi;

    double dSig_dr  = 2.0 * r;
    double dDel_dr  = 2.0 * r - 2.0 * M;
    double term_r = (dDel_dr*Sigma - Delta*dSig_dr)/(Sigma*Sigma) * s.pr*s.pr 
                  - dSig_dr/(Sigma*Sigma) * s.ptheta*s.ptheta;

    double num_gtt = (r2 + a2)*(r2 + a2);
    double d_num_gtt = 4.0*r*(r2+a2);
    double d_inv_gtt_dr = - ( (d_num_gtt*Delta - num_gtt*dDel_dr)/(Delta*Delta) * Sigma - (num_gtt/Delta - a2*sin2)*dSig_dr ) / (Sigma*Sigma);

    double num_gtp = a*(r2 + a2);
    double d_num_gtp = 2.0*a*r;
    double d_inv_gtp_dr = - ( (d_num_gtp*Delta - num_gtp*dDel_dr)/(Delta*Delta) * Sigma - (num_gtp/Delta - a)*dSig_dr ) / (Sigma*Sigma);

    double term_r_phi = 0.0;
    if (!zero_L) {
        double d_inv_gpp_dr = ( (- ( -a2*dDel_dr/(Delta*Delta) ) * Sigma) - (1.0/sin2 - a2/Delta)*dSig_dr ) / (Sigma*Sigma);
        term_r_phi = d_inv_gpp_dr * pphi2;
    }

    d.pr = -0.5 * (
        term_r +
        d_inv_gtt_dr*s.pt*s.pt +
        term_r_phi +
        2.0*d_inv_gtp_dr*s.pt*s.pphi
    );

    double dSig_dth = -2.0 * a2 * cost * sint;
    double term_th = -(Delta*dSig_dth)/(Sigma*Sigma) * s.pr*s.pr
                   - dSig_dth/(Sigma*Sigma) * s.ptheta*s.ptheta;

    double d_inv_gtt_dth = - ( (-a2*2.0*sint*cost)*Sigma - (num_gtt/Delta - a2*sin2)*dSig_dth ) / (Sigma*Sigma);
    double d_inv_gtp_dth = - ( - (num_gtp/Delta - a)*dSig_dth ) / (Sigma*Sigma);

    double term_th_phi = 0.0;
    if (!zero_L) {
        double d_csc2 = -2.0 * cost / (sin2 * sint); 
        double d_inv_gpp_dth = ( (d_csc2)*Sigma - (1.0/sin2 - a2/Delta)*dSig_dth ) / (Sigma*Sigma);
        term_th_phi = d_inv_gpp_dth * pphi2;
    }

    d.ptheta = -0.5 * (
        term_th +
        d_inv_gtt_dth*s.pt*s.pt +
        term_th_phi +
        2.0*d_inv_gtp_dth*s.pt*s.pphi
    );

    d.pt = 0.0;
    d.pphi = 0.0;
}

// ================ RK4 ==================
void rk4_step(State& s, double h) {
    State k1, k2, k3, k4, t;
    
    get_derivatives(s, k1);

    t = s; t.r += k1.r*h*0.5; t.theta += k1.theta*h*0.5;
    t.phi += k1.phi*h*0.5; t.t += k1.t*h*0.5;
    t.pr += k1.pr*h*0.5; t.ptheta += k1.ptheta*h*0.5;
    get_derivatives(t, k2);

    t = s; t.r += k2.r*h*0.5; t.theta += k2.theta*h*0.5;
    t.phi += k2.phi*h*0.5; t.t += k2.t*h*0.5;
    t.pr += k2.pr*h*0.5; t.ptheta += k2.ptheta*h*0.5;
    get_derivatives(t, k3);

    t = s; t.r += k3.r*h; t.theta += k3.theta*h;
    t.phi += k3.phi*h; t.t += k3.t*h;
    t.pr += k3.pr*h; t.ptheta += k3.ptheta*h;
    get_derivatives(t, k4);

    s.r += h*(k1.r + 2*k2.r + 2*k3.r + k4.r)/6.0;
    s.theta += h*(k1.theta + 2*k2.theta + 2*k3.theta + k4.theta)/6.0;
    s.phi += h*(k1.phi + 2*k2.phi + 2*k3.phi + k4.phi)/6.0;
    s.t += h*(k1.t + 2*k2.t + 2*k3.t + k4.t)/6.0;
    s.pr += h*(k1.pr + 2*k2.pr + 2*k3.pr + k4.pr)/6.0;
    s.ptheta += h*(k1.ptheta + 2*k2.ptheta + 2*k3.ptheta + k4.ptheta)/6.0;
}

// ================= Main =================
int main() {
    std::cout << "Starting Kerr Ray Tracer...\n";
    std::cout << "Resolution: " << W << " x " << H << "\n";
    
    std::vector<double> image_data(W * H * 2);

    #pragma omp parallel for schedule(dynamic, 10)
    for (int y = 0; y < H; y++) {
        
        if (y % 20 == 0) {
            #pragma omp critical
            std::cout << "\rRendering row " << y << " / " << H << std::flush;
        }

        for (int x = 0; x < W; x++) {

            double alpha = ((double)x/W - 0.5) * FOV;
            double beta  = ((double)y/H - 0.5) * FOV * ((double)H/W);
            
            State s{};
            s.r = R_CAM;
            s.theta = THETA_CAM;
            s.pt = -1.0;
            s.pphi = -alpha * std::sin(THETA_CAM);
            s.ptheta = -beta;

            double r = s.r; double th = s.theta;
            double r2=r*r; double a2=a*a;
            double Sigma = r2 + a2*std::cos(th)*std::cos(th);
            double Delta = r2 - 2*M*r + a2;
            double sin2 = std::sin(th)*std::sin(th);
            
            double inv_grr = Delta / Sigma;
            double inv_gth = 1.0 / Sigma;
            double inv_gtt = - ((r2 + a2)*(r2 + a2)/Delta - a2*sin2) / Sigma;
            double inv_gpp = (1.0/sin2 - a2/Delta) / Sigma;
            double inv_gtp = - (a*(r2 + a2)/Delta - a) / Sigma;

            double V_angular = 
                inv_gtt*s.pt*s.pt + 
                inv_gpp*s.pphi*s.pphi + 
                2.0*inv_gtp*s.pt*s.pphi + 
                inv_gth*s.ptheta*s.ptheta;

            int escaped_status = 0; 
            double sky_brightness = 0.0;

            if (V_angular < 0.0) {
                s.pr = -std::sqrt(-V_angular / inv_grr);

                for (int i = 0; i < MAX_STEPS; i++) {
                    double dist_to_horizon = s.r - R_HORIZON;
                    
                    // Using larger steps when far to save performance but extremely small steps near the horizon for precision
                    double h;
                    
                    if (dist_to_horizon > 100.0) {
                        h = 10.0;
                    } else if (dist_to_horizon > 10.0) {
                        h = 0.5;
                    } else {
                        h = 0.01 * dist_to_horizon;
                    }
                    
                    // Pole Safety
                    double sin_theta = std::abs(std::sin(s.theta));
                    if (sin_theta < 0.1) h *= (sin_theta * 5.0); // Smooth scaling

                    // Strict clamps
                    if (h > 5.0) h = 5.0;
                    if (h < 1e-5) h = 1e-5;

                    rk4_step(s, h);

                    // Checks
                    if (s.r < R_HORIZON + 0.005) {
                        escaped_status = 0; // Shadow
                        break; 
                    }
                    if (s.r > MAX_DIST) { 
                        escaped_status = 1; // Sky
                        break; 
                    }
                }

                if (escaped_status == 1) {
                    // ===== SMOOTH SKY SHADER =====
                    
                    // sin^2 to keep it always positive and smooth
                    double stripes_v = std::sin(30.0 * s.phi);
                    double stripes_h = std::sin(3.0 * s.theta);
                    
                    // soft grid pattern
                    double grid = 0.5 + 0.5 * (stripes_v * stripes_h);
                    
                    // Sharpen slightly using power, but keep smooth
                    sky_brightness = std::pow(grid, 0.8);

                    // Doppler effect
                    double beaming = 1.0 - 0.5 * (alpha / (FOV/2.0)); 
                    sky_brightness *= beaming;

                    sky_brightness = 0.05 + 0.95 * sky_brightness;
                    
                    if (sky_brightness < 0.0) sky_brightness = 0.0;
                    if (sky_brightness > 1.0) sky_brightness = 1.0;
                }
            }

            image_data[2*(y*W+x)+0] = (double)escaped_status;
            image_data[2*(y*W+x)+1] = sky_brightness;
        }
    }

    std::cout << "\nWriting output to RToutput.dat...\n";
    std::ofstream out("RToutput.dat");
    out << W << " " << H << "\n";
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            out << x << " " << y << " " 
                << image_data[2*(y*W+x)+0] << " " 
                << image_data[2*(y*W+x)+1] << "\n";
        }
    }
    out.close();

    std::cout << "Done.\n";
    return 0;

}
