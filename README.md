# Digital markers of circadian disruption

This code computes digital markers of circadian disruption based on activity, sleep, and heart rate data from the wearable devices. These digital markers can be powerful tools for future endeavors in studying the impact of circadian and sleep disruption on mental health and depression risks. 

## Basics

The following algorithm computes three different digital markers of circadian disruption, each represent different aspects of the underlying circadian and sleep physiology:

- The misalignment between Circadian Rhythms in the Central Oscillator (i.e., CRCO) and the sleep-wake cycle. We name this as the `CRCO-sleep misalignment`.
- The misalignment between Circadian Rhythms in the Peripheral Oscillator (i.e., CRPO) and the sleep-wake cycle. We name this as the `CRCO-sleep misalignment`. We particularly pay close attention to the peripheral circadian rhythms in the cardian pacemaker, following our previous works analyzing CRPO in the heart rate [(Bowman et al., 2021)](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(21)00106-5), [(Huang et al., 2021)](https://www.frontiersin.org/journals/digital-health/articles/10.3389/fdgth.2021.727504/full), [(Mayer et al., 2022)](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(22)00118-5), and [(Kim et al., 2023a)](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2023.0030).
- The misalignment between two different internal circadian rhythms (i.e., the misalignment between CRCO and CRPO). We name this as the `Internal misalignment`.

Each of these digital markers have been shown to be negatively associated with the daily mood scores. In particular, the `CRCO-sleep misalignment` furthermore predicts the long-term depression risks among young, working populations measured by the Patient Health Questionairre 9 (PHQ-9). For more details on how each markers were computed, please see `Methods` section of [(Lee et al., 2024; TO BE PUBLISHED)](https://www.nature.com/npjdigitalmed/).

## Data

- Data used in the study were collected as part of the [Intern Health Study](https://www.internhealthstudy.org/). The de-identified data from Intern Health Study that supports our findings may be available from the corresponding author upon reasonable request.

- For more information about the data, please refer to [(Sen et al., 2010)](https://jamanetwork.com/journals/jamapsychiatry/fullarticle/210823) or [(Guille et al., 2015)](https://jamanetwork.com/journals/jamapsychiatry/fullarticle/2467822).

## Pre-processing the data

Users will first need to convert their raw heart rate, sleep, and steps data into three .csv files that can be used as an input to the algorithms, with the following specific format:

- Column 1 should contain timestamp information in which the data was collected, in Unix timestamp (i.e., seconds since standard epoch of 1/1/1970).
- Column 2 should contain data value (e.g., measured heart rate in bpm, step counts, or sleep annotations).
- Column 3 should contain information about the data source (e.g., Apple Watch, iPhone, or FitBit), if applicable.

Additionally, if available, users will need to create a separate `subject.txt` file containing informations about

- Subject ID
- Local time zone information (e.g., America/Detroit or America/Chicago).

Note that the codes can be executed even without these information. In such case, time information will be converted into UTC time, and the algorithm will be executed in UTC timestamp. We have provided an examplary, deidentified wearable data collected from the Apple Watch collected as part of our previous study [(Bowman et al., 2021)](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(21)00106-5).

## Notes

- Running the comprehensive analysis with the provided sample data should take only about a minute.

- Users will need to appropriately change `line 5`, `line 40`, `line 49`, and `line 60` in `main.m` with their appropriate local directory to successfully import preprocessed data.

- As part of the algorithm, four separate `results.csv` files will be generated:
    1. `Sleep_results.csv` containing the estimated timing of sleep onset, sleep midpoint, and sleep offset for each day. This is obtained by running the funtion `run_sleep_analysis.m`.
    2. `CRPO_results.csv` containing the estimated parameter values describing CRPO in heart rate for each day. This is obtained by running the funtion `ALSM_nonlinear_version.m` based on our pervious work [(Kim et al., 2023a)](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2023.0030).
    3. `CRCO_results.csv` containing the estimated parameter values describing CRCO for each day. This is obtained by running the funtion `LSKF_Circadian.m` based on our pervious work [(Kim et al., 2023b)](https://epubs.siam.org/doi/abs/10.1137/22M1509680).
    4. `Digital_circadian_disruption_markers.csv` containing the estimated daily misalignment values for each of the three markers described as above.

- Each `CRCO_results.csv` and `CRPO_results.csv` should contain additional information about the estimated CRCO and CRPO, respectively, such as the daily phase estimates, uncertainty in phase estimates, circadian amplitude, mesors, noise level, etc.

- `Digital_circadian_disruption_markers.csv` should contain four columns, corresponding to `Datetime (dd-mm-yy)`, `CRCO-sleep misalignment`, `CRPO-sleep misalignment`, and `Internal misalignment`. 
 
## Important Note

Please note that the circadian disruption markers (i.e., CRCO-sleep and CRPO-sleep misalignments) can be executed without necessarily having a minute-by-minute sleep data for running `run_sleep_analysis.m`. In such case, however, users can manually provide information about `sleep onset`, `sleep midpoint`, and `sleep offset` information collected from the device to calculate these markers. 

## License

This software is open source and under an MIT license.
