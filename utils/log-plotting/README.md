# Monitoring log files

This ICON utility provides visual monitoring of job performances and wind speed information retrieved from the job log file.

## In this README

- [Features](#features)
- [Usage](#usage)
  - [Using mkexp](#using-mkexp)
  - [Using single modules](#using-single-modules)
- [Custom modules](#custom-modules)
  - [Analyzing modules](#analyzing-modules)
  - [Plotting modules](#plotting-modules)

## Features

### Extraction of information from the ICON log files.
By default, three types of information are extracted from the log files and saved to csv-tables.
- The times of each computing process
- Wind speed information
- Cycle times

### Creating plots from the extracted data.
- The total times of each compute process are displayed by job and rank.
- Wind speed distributions and time evolutions
- The cycle times are plotted as simulated days per day (SDPD) of one experiment.
- The format of the plots can be configured in the [config file](run/standard_experiments/DEFAULT.config) under the section `[[mon_log]] plot_format:`.

### Displaying the figures on an HTML page
All created plots are displayed on an overview page. This `index.html` page is by default saved to `$MON_DIR/log_plotting`. On top of this overview page, there appear plots that have been prioritized in the configuration. These plots on the top of the page can be set in the [config file](run/standard_experiments/DEFAULT.config) under the section `[[mon_log]] prioritized_plots:`.

## Usage

Using log-monitoring on existing data:

### Using mkexp

1. Setup ICON runtime environment
```shell
EXP_ID=abc1234
PROJECT=ab1234
ICON_DIR=${HOME}/icon-mpim

WORK_DIR=/work/${PROJECT}/${USER}/$(basename ${ICON_DIR})
SCRIPT_DIR=${ICON_DIR}/experiments/${EXP_ID}/scripts
MON_DIR=${WORK_DIR}/monitoring/${EXP_ID}

cd ~/icon-mpim/build
module use ${ICON_DIR}/etc/Modules
module rm cdo
module add icon-levante
```
2. Select the date in the simulation according to the desired job ID
You can browse for the date or job ID in `${SCRIPT_DIR}/${EXP_ID}.log`.

3. Execute log monitoring
```shell
DATE=2020-01-01
cd ${SCRIPT_DIR}
editexp ${DATE}
./update
./${EXP_ID}.mon_log -c ${SCRIPT_DIR}/${EXP_ID}.dump ${DATE}
```
4. Open overview page `${MON_DIR}/log_plotting/index.html`

### Using single modules

To extract or plot data from the regular **job log file** use [`run_log_parser.py`](utils/log_monitoring/run_log_parser.py) or [`run_log_plotter.py`](utils/log_monitoring/run_log_plotter.py) and use the option `--custom_modules` to add the modules which you want to get executed.
```shell
python run_log_parser.py ${SCRIPT_DIR}/${EXP_ID}.run.${JOB_ID}.log -o <output directory> --job_id ${JOB_ID} --exp_id ${EXP_ID} --custom_modules AnalyzeTimer AnalyzeWind
```
```shell
python run_log_plotter.py <csv-file directory> -o <output directory> --custom_modules PlotTimer PlotWindSpeeds --plot_format <plot_format>
```

For the **experiment status log** use [`exp_status_log_parser.py`](utils/log_monitoring/exp_status_log_parser.py) to extract and [`PlotSDPD.py`](utils/log_monitoring/PlotSDPD.py) to plot the cycle times.
```shell
python exp_status_log_parser.py ${SCRIPT_DIR}/${EXP_ID}.log -o <output directory> -
```
```shell
python PlotSDPD.py <csv-file directory> -o <output directory> --plot_format <plot_format>
```

To generate the **overview page** use [`build_index_html.py`](utils/log_monitoring/build_index_html.py) with the directory containing all plots you want to display.
```shell
python build_index_html.py <plot directory> --exp_id ${EXP_ID} --plot_format <plot_format>
```

## Custom modules

For better customisability, the possibility of adding your own modules is provided. These can be either classes to extract data from the log file or to create figures out of data from existing extracted data tables. Follow the logic below to integrate custom modules into log monitoring.

### Analyzing modules
In case you would like to extend the log monitoring on data in the log file you are interested in, these are the steps to follow:

1. Write a custom module:

To work with the `run_log_parser.py` there are a few requirements. The superclass [`BaseAnalyzer`](utils/log_monitoring/BaseAnalyzer.py) provides the general structure of each analyzer class and some default functions. Import and use this superclass when writing your analyzer class. It will raise an error if you forget to implement the two main functions `analyze_line` and `processing`. The function `analyze_line` is supposed to analyze the data of interest from a given line of the log file. In the next step `processing` creates a data frame from the collected data. Please note that conventionally, the last two columns should be the experiment name and job ID. The plotting modules rely on this convention to work without the need to supply experiment and job ID.

In principle, your custom module should be structured like this:

```python
class AnalyzeSomething:
    def __init__(self):
        super().__init__()
        self.data = []
        self.filename = 'Something'

    def analyze_line(self, line):
        #extract the data of interest
        return self.data

    def processing(self, exp_id, job_id):
        #create a pandas dataframe out of the collected data
        self.df = pd.DataFrame(self.data, columns=("Column1", "Column2"))
        #last columns should be the experiment name and job ID
        self.df["exp_name"] = [exp_id] * len(self.df)
        self.df["job_id"] = [job_id] * len(self.df)
        return self.df

custom_analyzers = [AnalyzeSomething(), ]
```
It is important to name the functions exactly `analyze_line` and `processing` as well as giving a `filename` under which the csv table is saved in the end.

Your python file should also include a variable called `custom_analyzers`. It is supposed to be a list of the analyzers defined in your python file. This step is necessary so the new analyzing class gets recognized by `run_log_parser.py`.

2. Set up config:

Go to [`/icon-mpim/run/standard_experiments/DEFAULT.config`](run/standard_experiments/DEFAULT.config) and add the name of your module (without `.py`ending) to `custom_modules` under `[[mon_log]]`.

This step completes the setup of a new analyzing module. The newly generated csv-files should show up together with the other output in the same directory after the first execution.

### Plotting modules

Adding custom plot modules follows the same approach as for the analyzers. The superclass [`BasePlotter`](utils/log_monitoring/BasePlotter.py) provides the general structure and ready to use function called `extract_data`. If the function `plot` is not defined in the new plotter subclass a NotImplementedError will be raised. The `plot` function must include saving the plots.

As for the analyzers, a file containing the plotter should also include the variable `custom_plotter`.
