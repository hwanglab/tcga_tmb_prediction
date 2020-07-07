First, please go to website: https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%220f4fb205-b13c-4a85-97b3-779429e6ccdd%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22Diagnostic%20Slide%22%5D%7D%7D%5D%7D&searchTableTab=files

Then, download the example pathlogy slide: TCGA-BL-A3JM-01Z-00-DX1.33E53972-CEA4-4D84-A5D2-7DAD7B0C27F8.svs
we cannot upload any slide in github due to size limiatation.


This is a simple tmb_prediction example for a small whole slide image.

You can get the idea how we make tmb predictions in our pipeline.

For more detail experiments, you need to go through our github steps.

This example has 3 steps.
step1: select representative tumor tiles (matlab code)
step2: extract image features (python code)
step3: tmb prediction (matlab code)

You can run step1-step2-step3
Or if you only have matlab:
You can also run step3 directly

For questions: mxu@ualberta.ca
