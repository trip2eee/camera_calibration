# Camera Calibration
Camera calibration tool.

## Axis
### Image coordinate
* $u$: Lateral axis. Leftward direction is positive.
* $v$: Vertical axis. Downward direction is positive.

### World coordinate
* $x$: Lateral axis. Leftward direction is positive
* $y$: Vertical axis. Downward direction is positive.
* $z$: Longitudinal axis. Forward direction is positive.


## Intrinsic camera parameters

* $f_u$ focal length in u-axis
* $f_v$ focal length in v-axis
* $cu$ center of image in u-axis
* $cv$ center of image in v-axis
* skew is assumed to be 0

$$K = \begin{bmatrix} f_u & 0 & cu \\\\ 0 & f_v & cv \\\\ 0 & 0 & 1 \end{bmatrix}$$

## Extrinsic camera parameters
* $t_x$ translation in x-axis
* $t_y$ translation in y-axis
* $t_z$ translation in z-axis

* $\phi$ rotation around x-axis
* $\theta$ rotation around y-axis
* $\psi$ rotation around z-axis

## Lens distortion

This code computes radial distortion only. Tangential distortion is ignored.

Lens distortion model converts undistorted image coordinate $(u,v)$ into distorted image coordinate $(u_d, v_d)$


$u_n = \frac{(u - cu)}{f_u}$

$v_n = \frac{(v - cv)}{f_v}$

$r_n^2 = u_n^2 + v_n^2$

$ratio = \kappa_1 r_n^2 + \kappa_2 r_n^4 + + \kappa_3 r_n^6$

$u_d = f_u \cdot ratio \cdot u_n + cu$

$v_d = f_v \cdot ratio \cdot v_n + cv$


## Initial value

Compute homography with $z=0$

$$\begin{bmatrix}w\cdot u\\\\w\cdot v\\\\w\end{bmatrix} = H\begin{bmatrix}x \\\\ y \\\\ 1\end{bmatrix}$$

where $H = \lambda K\begin{bmatrix}r_1 & r_2 & t\end{bmatrix}$, $r_1$ and $r_2$ are the 1st and the 2nd column vectors of the rotation matrix $R$ and $t$ is translation vector. $\lambda$ is an arbitrary scalar.


Let's denote $H = \begin{bmatrix}h_1 & h_2 & h_3\end{bmatrix}$. Then $K$, $r_1$ and $r_2$ can be decomposed by exploiting orthogonality of $R$.

$h_1=\lambda K r_1$

$h_2=\lambda K r_2$

$r_1 = \frac{1}{\lambda} K^{-1}h_1$

$r_2 = \frac{1}{\lambda} K^{-1}h_2$

$r_1^T r_2 = h_1^TK^{-T}K^{-1}h_2 = 0$

$r_1^T r_1 = r_2^T r_2 = h_1^TK^{-T}K^{-1}h_1 = h_2^TK^{-T}K^{-1}h_2$

The above equations can be rewritten as

$h_1^TBh_2 = 0$

$h_1^TBh_1 - h_2^TBh_2 = 0$

$$ B = K^{-T} K^{-1} = \begin{bmatrix}\frac{1}{f_{u}^{2}} & 0 & - \frac{cu}{f_{u}^{2}}\\\\0 & \frac{1}{f_{v}^{2}} & - \frac{cv}{f_{v}^{2}}\\\\- \frac{cu}{f_{u}^{2}} & - \frac{cv}{f_{v}^{2}} & \frac{cu^{2}}{f_{u}^{2}} + \frac{cv^{2}}{f_{v}^{2}} + 1\end{bmatrix} = \begin{bmatrix}b_{11} & 0 & b_{13} \\\\ 0 & b_{22} & b_{23} \\\\ b_{13} & b_{23} & b_{33} \end{bmatrix}$$

Then the elements of $B$ matrix can be computed by the following homogeneous equation.
$$Vb = 0$$
where
$$V=\begin{bmatrix}h_{11} h_{12} & h_{11} h_{32} + h_{12} h_{31} & h_{21} h_{22} & h_{21} h_{32} + h_{22} h_{31} & h_{31} h_{32}\\\\h_{11}^{2} - h_{12}^{2} & 2 h_{11} h_{31} - 2 h_{12} h_{32} & h_{21}^{2} - h_{22}^{2} & 2 h_{21} h_{31} - 2 h_{22} h_{32} & h_{31}^{2} - h_{32}^{2}\end{bmatrix}$$
$$b = \begin{bmatrix}b_{11} & b_{13} & b_{22} & b_{23} & b_{33}\end{bmatrix}^T$$

For detailed process of symbolic computation, please refer to [initial_value.ipynb](initial_value.ipynb)

# Reference
Z. Zhang. A flexible new technique for camera calibration. IEEE Transactions on Pattern Analysis and Machine Intelligence, 22(11):1330-1334, 2000

