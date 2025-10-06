# Generation-of-Samples-for-Target-Years-Using-Ensemble-Historical-Machine-Learning-Models
Generation of Samples for Target Years  Using Ensemble Historical Machine Learning Models




Due to the non-public nature of Nenjiang data, the research code in this study uses the Indiana research area in the United States as an example 

First:
Run 2020_sample_generation.js 2021_sample_generation.js 2022_sample_generation.js
and acquire historical samples:
maizesample_purify_2020_withfeature,maizesample_purify_2021_withfeature,maizesample_purify_2022_withfeature,soybeansample_purify_2020_withfeature,soybeansample_purify_2021_withfeature,soybeansample_purify_2022_withfeature,othercropsample_purify_2020_withfeature,othercropsample_purify_2021_withfeature,othercropsample_purify_2022_withfeature.

Second:
Input the samples generated in the previous step into file:transfer_model_and_generate_target_year_sample.js
and acquire target year samples:
maizesample_2023,soybeansample_2023,othercropsample_2023

Third:
Input the target year samples generated in the previous step into file:classification_in_2023.js
and classification.
