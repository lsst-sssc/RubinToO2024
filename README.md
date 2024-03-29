# RubinToO2024
Rubin ToO 2024: Envisioning the Vera C. Rubin Observatory LSST Target of Opportunity program for Solar System Science

This is a fork of https://lssttooworkshop.github.io with a small amount of Solar System-specific added features.

The notebook plots_RubinToO2024.ipynb has the following features:
- visualization of the ToO observing strategies
- light curves plots
- exposure time calculator (ETC) to convert exposure times to limiting magnitudes and vice-versa and modified to account for the near-Sun twilight survey and trailing losses
- time budget calculations


## to use on Colab:
we provided a Google Colab notebook version of this notebook at https://drive.google.com/file/d/10fRrqK9eDnaR2Va7RoYAkSlaVEAr5zMb/view?usp=sharing


To use this notebook: 
- open the link 
- click on "Open with Colab" (if that is not your default at this point you will see the ugly json file)
- save a copy of the notebook in your Drive (you do not have permission to edit the notebook directly)
<img width="446" alt="Screen Shot 2024-03-18 at 11 23 50 AM" src="https://github.com/fedhere/RubinToO2024/assets/1696902/da4af5ff-5d9f-4d1b-92c3-32f09cf00dd5">

- run the top cells of the notebook to connect mount your google drive, you will have some pop-up windows to go through at this point and give permissions
<img width="272" alt="Screen Shot 2024-03-18 at 11 14 09 AM" src="https://github.com/igorandreoni/RubinToO2024/assets/1696902/32788762-cd39-42d8-b3b3-11f4adc16234">

- run the cells that clone the original repo

- you will only need to run this once, after the first time you run the notebook 
<img width="629" alt="Screen Shot 2024-03-18 at 11 25 53 AM" src="https://github.com/fedhere/RubinToO2024/assets/1696902/79b8d46d-acd8-4ec2-9352-5e7671eae572">

## running tests
since this is not fully packaged, you will need to run the tests as:
```
python -m pytest
```
in order to put the current directory onto the `PYTHONPATH` and for the tests to find the code.
