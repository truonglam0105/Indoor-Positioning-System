Description
This project develops an indoor positioning system (IPS) using Wi-Fi signal strength (RSSI) to predict device location inside a building.

Dataset

6,600+ Wi-Fi signal measurements collected across 166 indoor locations

Multiple device orientations for robustness testing

Methods

K-Nearest Neighbors (KNN) Regression

Trilateration with Least Squares Regression

Noise filtering and preprocessing for signal stability

Results

KNN achieved 3.18m RMSE, outperforming Least Squares (4.81m RMSE)

Accuracy influenced by signal attenuation, device orientation, and obstacles

How to Run

Implemented in Python (scikit-learn, numpy, pandas)

Run main.py with training dataset to reproduce results

Future Work

Integrate deep learning methods (CNN/RNN) for improved signal modeling

Extend to multi-floor indoor positioning
