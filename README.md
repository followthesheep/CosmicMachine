# Investigating cosmic ray and bad pixel corrections

## Files

- reduced/s170517_a018002_Kn3_035_nomask.fits - ran with no mask

## Results

We parameterize the classification according to the ``classification_report`` of scikit-learn. The parameters are:

The precision is the ratio tp / (tp + fp) where tp is the number of true positives and fp the number of false positives. The precision is intuitively the ability of the classifier not to label as positive a sample that is negative.

The recall is the ratio tp / (tp + fn) where tp is the number of true positives and fn the number of false negatives. The recall is intuitively the ability of the classifier to find all the positive samples.

The F-beta score can be interpreted as a weighted harmonic mean of the precision and recall, where an F-beta score reaches its best value at 1 and worst score at 0.

The F-beta score weights recall more than precision by a factor of beta. beta == 1.0 means recall and precision are equally important.

The support is the number of occurrences of each class in y_true


- Split the pixels into two categories: cosmic rays, valid pixels

```
precision    recall  f1-score   support

Cosmic Ray       0.99      0.88      0.93       500
Valid Data       0.94      0.99      0.97      1000

avg / total       0.96      0.96      0.96      1500

```
This is much better at distinguishing the valid pixels from the cosmic ray pixels. There are however about 12% of cosmic ray pixels that are marked as good (recall) and about 1% of valid data is marked as cosmic rays.

- Split the pixels into three groups: cosmic rays, valid observed pixels, and OH line pixels

```
                  precision    recall  f1-score   support

 Cosmic Ray       0.94      0.93      0.94       500
 Valid Data       0.49      0.96      0.65       500
   OH lines       0.33      0.01      0.02       500

avg / total       0.59      0.63      0.53      1500
```

This has trouble with identifying OH lines from regular data points
