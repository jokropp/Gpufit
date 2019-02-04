#ifndef GPUFIT_PLISIN_CUH_INCLUDED
#define GPUFIT_PLISIN_CUH_INCLUDED
#include <cmath>
/*
* 2 Dimensional Fit of PLI Intensities function
* Paramters are
    * phi - direction
    * alpha - inclination
    * t_rel - t_rel
* Variables are
    * rot - rotation
    * tilt - tilting direction (None, 0, 90, 180, 270) -> (0, 1, 2, 3, 4)
* User Info needs to be
[tau, sub_sample, n_tilts, n_rots, tilt_1, ..., tilt_n, rot_1, ..., rot_n]
*/
const REAL pi = 3.14159265358979f;

__device__ void calculate_plisin(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    // extract tau, sub_sample, rotation angle and tilt direction
    REAL * user_info_float = (REAL*) user_info;
    REAL tau = user_info_float[0];
    int sub_sample = int(user_info_float[1]);
    int n_tilts = int(user_info_float[2]);
    int n_rots = int(user_info_float[3]);
    // apply sub sampling
    int index = point_index;
    index *= sub_sample;
    REAL tilt = user_info_float[(int) index/n_rots + 4];
    REAL rot = user_info_float[index%n_rots + n_tilts + 4];

    // get params (tilt them if needed)
    REAL phi, alpha, t_rel;
    phi = parameters[0];
    alpha = parameters[1];
    t_rel = parameters[2];

    if (tilt == 0.0f){

    }else{
        REAL angles[4] = {0.0f, pi/2, pi, 3*pi/2};
        REAL psi_t = angles[int(tilt) - 1];

        // create rotation_matrix
        REAL R[3][3] = {0};
        R[0][0] = cos(tau) * cos(psi_t)*cos(psi_t) + sin(psi_t)*sin(psi_t);
        R[0][1] = (cos(tau) - 1) * sin(psi_t) * cos(psi_t);
        R[0][2] = cos(psi_t) * sin(tau);

        R[1][0] = (cos(tau) - 1) * sin(psi_t) * cos(psi_t);
        R[1][1] = cos(tau) * sin(psi_t)*sin(psi_t) + cos(psi_t)*cos(psi_t);
        R[1][2] = sin(psi_t) * sin(tau);

        R[2][0] = -cos(psi_t) * sin(tau);
        R[2][1] = -sin(psi_t) * sin(tau);
        R[2][2] = cos(tau);

        // create orientation vector
        REAL vec[3] = {0};
        vec[0] = cos(alpha) * cos(phi);
        vec[1] = cos(alpha) * sin(phi);
        vec[2] = sin(alpha);

        // multiply vector with rotation matrix
        REAL temp_vec[3] = {vec[0], vec[1], vec[2]};
        for (int i = 0; i < 3; i++){
            REAL sum = 0.0f;
            for (int j = 0; j < 3; j++){
                sum += temp_vec[j] * R[i][j];
            }
            vec[i] = sum;
        }

        // get params from orientation vector
        phi = atan(vec[1] / REAL(vec[0]));
        alpha = asin(vec[2]);
        t_rel /= cos(tau);
    }

    // value
    value[point_index] = sin(2 * (rot - phi)) * sin(pi/2 * t_rel * cos(alpha) * cos(alpha));

    // derivatives
    REAL * current_derivatives = derivative + point_index;
    REAL cosa2 = cos(alpha) * cos(alpha);
    REAL ret = pi/2 * t_rel * cosa2;
    REAL temp = sin(2*(rot - phi));
    current_derivatives[0 * n_points] = -2 * sin(ret) * cos(2 * (rot - phi)); //phi
    current_derivatives[1 * n_points] = -temp * cos(ret) * pi * t_rel * cos(alpha) * sin(alpha); // alpha
    current_derivatives[2 * n_points] = temp * pi/2 * cosa2 * cos(ret); //t_rel
}

#endif
