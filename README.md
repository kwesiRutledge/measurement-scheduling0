# measurement-scheduling0
This repository contains the joint work of Antoine Aspeel and Kwesi Rutledge on a self-triggering control problem. Their work part of a joint project between Prof. Ozay (University of Michigan) and Profs. Macq and Jungers (Université catholique de Louvain).

## How this Repository is Organized

In `examples`, the code which is used to produce the examples in this repositories associated paper is written.

In `lib`, the julia helper functions which allow the examples to run is written.

In `tests`, the testing suite for the helper functions is contained.

## Problem Statement
This project aims to control systems of the form
$$
\begin{array}{rll}
    x_{t+1} & = A x_t + B u_t + k + w_t, & w_t \in \mathcal{W} \\ \cr
	y_t & =
	    \begin{cases}
		    C x_t + v_t, & \sigma^m_t = 1 \\ \cr
			\emptyset, & \sigma^m_t = 0
		\end{cases}
		& v_t \in \mathcal{V} \\ \cr
	z_t &=D x_t+d 
\end{array}
$$

(rest of problem will be added later)