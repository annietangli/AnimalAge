# AnimalAge
1. code: report.R
2. rmarkdown: report.Rmd
3. pdf: report.pdf

Strongly recommend to ***download the pdf*** (instead of viewing on GitHub) to see plots rendered correctly in high resolution

## A short summary
In this machine learning project, I'll solve a regression problem of predicting animal max-longevity, which is the recorded highest lifespan for each species (e.g. Bowhead whale: 211 years). Since this is a regression problem, the goal is to achieve a low rmse (root mean squared error). The requirement is to develop at least 2 models with 1 model more advanced than linear regression.

My approach is to create ensembles: combining results of different methods into one that hopefully improves the result. I developed two ensemble models. The first ensemble is linear, combining cubist and lm. The second ensemble is non-linear, combining rf and rpart. The results of both ensembles are similar: linear achieved a rmse of 8.31 years, and non-linear achieved a rmse of 8.91 years. I think my results are good considering animals max-longevity could range from about 2 to 200 years.
