# Unscented Kalman Filter Project 
Self-Driving Car Engineer Nanodegree Program

---
[ukf_sample_1]: ./doc/sample-1-std_a_1.2_std_yawdd_0.5.png
[ukf_sample_2]: ./doc/sample-2.png
[ukf_sample_3]: ./doc/ukf-obj-pose.png


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Experiments result
![alt_text][ukf_sample_1]
![alt_text][ukf_sample_2]
![alt_text][ukf_sample_3]

## Generating Additional Data

see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

