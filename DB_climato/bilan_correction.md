correction with absorption
================
Flavien
15/11/2019

# The observation

We can look at match up data between the fluorescence of the first
BGC-Argo float profile and HPLC measurements. This insight indicate the
bias we have when we estimate Chla concentration from fluorescence.
Thus, the ratio \[Chla\]<sub>Fluo</sub>/\[Chla\]<sub>HPLC</sub> give us
the factor by which we have to divide the fluorescence estimation to fit
with HPLC Data. Currently the factor we use is the same for the whole
ocean and is not taking time or depth in consideration (Roesler *et al.*
2017). <br> So we decided to present a new way of correcting the
estimation of \[Chla\] by fluorescence that could take the variability
of the bias into account. <br> From HPLC data we can compute the
absorbtion of photosynthetical pigments a<sub>PS</sub>, by multiplying
their concentration by their specific absorbtion in solution (Data from
Bricaud and Claustre). Since the fluorometer excite at 470nm and the
Chla maximum absorbtion is at 440nm we assume that the ratio between
a<sub>PS</sub>440/a<sub>PS</sub>470 could be an index of the bias.<br>
![](Figs/unnamed-chunk-2-1.png)<!-- -->

<br> Here we see that the a<sub>PS</sub>440/a<sub>PS</sub>470 ratio is
effictively linked to the bias of the estimation of \[Chla\] by
fluorescence. Though, it can be a bias of visualisation because we group
data by bioregion and we mean it.<br> We can look at it by a more
precise way : ![](Figs/unnamed-chunk-3-1.png)<!-- --> <br> We show here
a non linear correlation between the two variables. This mean that we
can estimate the bias (aka the correction factor) from the ratio of
photosynthetical absorbtion.<br> From there we try to correct the
fluorescence, which give us this result :
![](Figs/unnamed-chunk-4-1.png)<!-- --> <br> This result leads us to
think that it is possible to rely on the
a<sub>PS</sub>440/a<sub>PS</sub>470 ratio to correct the \[Chla\]
estimation from fluorescence. <br> We need a way to estimate this ratio
anywhere the float isconsidering time and depth. <br> For this we choose
the classification method called Support Vector Machine. This tool will
be trained with a global HPLC database that is coupled with
climatological data. The choice of climatological data has been made to
be sure to conserv all our HPLC sample and not loose any information due
to cloud coverage. <br> The map of HPLC observation is the one below.

    ## [1] 1719    5

![](Figs/unnamed-chunk-5-1.png)<!-- -->

<br> The Data we included are listed below :<br>

<li>

rrs667

</li>

<li>

rrs555

</li>

<li>

rrs488

</li>

<li>

rrs443

</li>

<li>

rrs412

</li>

<li>

MLD

</li>

<li>

PAR

</li>

<li>

Depth

</li>

<br> Those are climatological values that have been computed every month
on a grid of 9km (downloaded
<a href="https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Monthly_Climatology/9km/">Here</a>)<br>
A quality controle have been applied on HPLC data :

<li>

Only lov HPLC data have been chose

</li>

<li>

Profiles with more than 4 points

</li>

<li>

Profiles with the first measurement in the first 10m and the last one
below the Zeu

</li>

<li>

Profiles with more than 500m depth

</li>

<li>

Profiles that are not in the argo/HPLC matchup databse

</li>

<br> Finally we have 1728 profiles in our databse, including 1382
profiles for training set and 346 for the validation.<br>

# SVM model

![](Figs/unnamed-chunk-6-1.png)<!-- -->

So the prediction by the SVM is good with a nice fit and a relatively
small error (RMSE). Nevertheless we can observe some outlayers that are
not well predicted by our model. Anyway we will now apply this model on
our database of matchup between HPLC measurement and Argo float
fluorescence. Then, from this prediction of the
a<sub>PS</sub>440/a<sub>PS</sub>470 ratio we will try to correct the
\[Chla\] estimation with the formula we fitted previously.<br>
![](Figs/unnamed-chunk-7-1.png)<!-- --> We observe there a low R² and a
significativ error, certainly due to outlyers that we can spot on the
plot. Nevertheless we can see that a part of the values have been quite
well predicted. Our quantification of error is sensitiv to outliers.<br>
We can still test a correction of the fluorescence based on this
prediction to see what happen. <br> \# Correction
![](Figs/unnamed-chunk-8-1.png)<!-- -->![](Figs/unnamed-chunk-8-2.png)<!-- -->

The RMSE of the corrected estimation of Chla is 0.26 whereas the RMSE of
non corrected estimation of Chla is 0.55. This means that we improved
the estimation of Chla. BUT do we really take the variabilty into
account ?<br> To test that we can first check the distribution of our
correction factor, and then test the quality of correction by a
constent.<br> ![](Figs/unnamed-chunk-9-1.png)<!-- -->

We have a distribution that looks normal and centered around 3. One
problem is that the minimum value of the correction factor is 1.6,
despite the fact that some region have an underestimation of Chla by
fluorescence. So this correction is not suitable for oligotrophic region
that have a Fluo/Chla close/or lower to 1… <br> The second problem can
be simply observed, by the performance of a correction by a factor 3,
which correspond to the center of the distribution of our normal
distribution.

![](Figs/unnamed-chunk-10-1.png)<!-- -->

The RMSE of the correction by this factor is 0.32, so it also improve
the quality of our estimation. This seems to indicate that our method
improve the \[Chla\] estimation because the correction factor we
calculate is between 1.6 and 5 which is reasonable for the
overestimation we are facing.<br> Finally, to picture it we can plot the
Fluo/\[Chla\] ratio by region before and after the correction.<br>
![](Figs/unnamed-chunk-11-1.png)<!-- --> <br>Most of the regions have a
Fluo/Chla value lower than 1 now, which show that our method use a
correction factor that is too high and not adapted to all regions.
