import os
import warnings
import re
import pickle
import time

import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from tqdm import tqdm
import uncertainties as unc
from uncertainties import unumpy as unp
# File containing useful redundant functions
import helpers

try:
    import tkinter as tk
    from tkinter import filedialog as fd
except ModuleNotFoundError:
    print("You do not have Tkinter installed, you cannot use the GUI :/")

class Reduced:
    """
    Object of Analysis of Neptune plus LA-MC-ICPMS data files in Boron analysis.
    """
    def __init__(self, rawfile=None, logfile=None):

        self.rawfile = rawfile
        self.logfile = logfile

    def rawfilefil(self):
        """Filter the raw files for analysis."""
        nfiles = False
        if len(self.rawfile) > 1 and type(self.rawfile) is tuple:
            # More than one data file (chosen with Tk)
            # Get all the file names
            names = [re.findall("([A-z]{1,})\d{1,}\.", f)[0] for f in self.rawfile]
            self.directory = os.path.dirname(self.rawfile[0])
            newfilename = self.directory + "/" + names[0] + "_r.tab"
            nfiles = True

        elif len(self.rawfile) == 1 and type(self.rawfile) is tuple:
            self.rawfile = self.rawfile[0]  # Extract the single filename from list
            self.directory = os.path.dirname(self.rawfile)
            newfilename = self.rawfile[:-4] + "_r.tab"
        else:
            self.directory = os.path.dirname(self.rawfile)
            newfilename = self.rawfile[:-4] + "_r.tab"

        ######################################################################
        # Check if the filtered file is already there
        ######################################################################
        if nfiles:  # If there are several files
            with open(newfilename, "w") as nf:
                # Loop through all the files and append content to newfile
                for i, rn in enumerate(self.rawfile):
                    with open(rn, "r") as rf:
                        _ = [rf.readline() for x in range(16)] # Skip header (16 lines on Neptune)
                        if i == 0:
                            nf.writelines(rf.readline())
                        else:
                            _ = rf.readline()
                        # All the lines starting with a number are for data.
                        # Filters the exo file to just have the data rows
                        nf.writelines(filter(helpers.float_check, rf.readlines()))

        else:  # For a single file.
            with open(self.rawfile, "r") as rf:
                # Opens the raw .exp file
                header = [rf.readline() for x in range(16)]  # Saves the header (Neptune specific)
                with open(newfilename, "w") as f:
                    # Write the column names (17th row)
                    f.writelines(rf.readline())
                    # Apply the filter to each line
                    f.writelines(filter(helpers.float_check, rf.readlines()))

        #######################################################################
        # Reading the data in and filtering unecesary columns
        #######################################################################
        self.lindata = pd.read_csv(newfilename, sep="\t")       
        ###### Filter the columns
        for i in self.lindata.columns:
            if "(" in i:
                # The Neptune software uses brackets to report automatic ratios
                self.lindata = self.lindata.drop([i], axis=1)
            elif "Unnamed" in i:
                # Empty columns are named Unnamed by pandas
                self.lindata = self.lindata.drop([i], axis=1)
        self.lindata = self.lindata.drop(["Cycle"], axis=1)

        # Remove the negative values
        self.lindata.iloc[:, 1:] = self.lindata.iloc[:, 1:].clip(lower=0)
        # Get rid of empty rows
        self.lindata = self.lindata.loc[~pd.isna(self.lindata["Time"])]
        # Save the data
        # self.lindata.to_csv(newfilename[:-3] + "csv", index=None)
        # Remove the tab delimited file
        os.remove(newfilename)


    def tkfile(self):
        """
        Method for having interactive Tk filedialog to choose the raw and log
        file instead of directly giving it.
        """
        tk.Tk().withdraw()  # Empty window of Tk
        self.rawfile = fd.askopenfilenames(title="Select the raw file(s)",
                                           filetypes=[("Raw files", "*.exp"),
                                                      ("All files", "*")])
        self.logfile = fd.askopenfilenames(title="Select the logfile(s)",
                                           filetypes=[("Log files", "*.csv")])

    def standard_detect(self):
        """Find the standards."""
        # Quite an ugly way of doing it, also not flexible 
        # TODO use regex to make a bit more flexible
        stdlist = ["nist612", "macs", "jct", "jcp", "blue", "yellow", "uwc-1", 
                   "uwc_1","uwc1", "uwc-3", "uwc_3", "uwc3", "kcsp"]
        
        self.stds = []
        for s in stdlist:
            l = [i for i, v in enumerate(self.sn) if re.findall(s, v)]
            if l:
                self.stds.append(s)
                        

    def spotsize(self):
        """Find the different spotsizes in the logfile."""
        # Make sure that there are no trailing white spaces
        self.rawlog.columns = [x.strip(" ") for x in self.rawlog.columns]
        # Look through the log file and correct the spot size if it is not in the
        # right format (i.e "40.0x40.0" instead of "40")
        for i in range(len(self.rawlog)):
            try: # to convert to integer
                int(self.rawlog["Spot Size (um)"][i])
            except ValueError: # Correct the wrong format
                self.rawlog.loc[i, "Spot Size (um)"] = re.findall("(\d{1,})",
                                                                  self.rawlog["Spot Size (um)"][i])[0]
        # Make a list of the different spot sizes using the "set" function
        self.spots = list(set(self.rawlog["Spot Size (um)"]))

    def logtime(self, check=False, manstart=False):
        """Convert the time in the file into real time and create data dictionary."""
        #######################################################################
        # If the rawfile and logfile were not directly provided, prompt user
        # to choose them.
        #######################################################################
        if self.rawfile is None and self.logfile is None:
            # If the exp and log files are not given during the instantiation
            self.tkfile()  # Select the files with GUI
            self.rawfilefil()  # Filter the rawfiles
            self.rawlog = pd.read_csv(self.logfile[0])  # Read the logfile
        else:  # If the files were provided OR the method has already been run
            self.rawfilefil()
            if type(self.logfile) == tuple:  # If the files were chosen with GUI
                self.rawlog = pd.read_csv(self.logfile[0])
            else:
                self.rawlog = pd.read_csv(self.logfile)

        self.spotsize()  # Get the different spotsizes
        print(f"\n\nThis session has the following spotsizes: {self.spots} um\n\n")
        
        # Transform the Timestamp in log file as Unix timestamp
        self.rawlog["Timestamp"] = pd.to_datetime(self.rawlog.Timestamp)
        self.rawlog["Time"] = pd.to_datetime(self.rawlog.Timestamp).values.astype(np.int64)

        #######################################################################
        # Go through the log file for sample names, start and end of laser.
        #######################################################################
        
        timela_start, timela_end, self.sn = [], [], []
        self.spotsizes = []
        self.fluences = []
        self.rep_rates = []
        self.xpos = []
        self.ypos = []
        # TODO cleanup this mess, maybe read log in "helpers" and return a DF
        # of all the parameters?
        # loop through the log file (ONLY WORKS FOR GEOSTAR software)
        for i in range(len(self.rawlog)):
            if not pd.isna(self.rawlog.iloc[i, 4]):  # If the cell is not empty
            
                # Loop through the log file to extract start and end time of each sample
                self.spotsizes.append(int(self.rawlog.iloc[i+3, 13]))
                self.sn.append(self.rawlog.iloc[i, 4].lower())
                self.fluences.append(self.rawlog.iloc[i+3, -2])
                self.rep_rates.append(self.rawlog.iloc[i+3, 11])
                self.xpos.append(self.rawlog.iloc[i+3, 5])
                self.ypos.append(self.rawlog.iloc[i+3, 6])
                timela_start.append(self.rawlog.iloc[i+3, -1])
                timela_end.append(self.rawlog.iloc[i+4, -1])

        data = self.lindata
        try:
            # If the timestamp is already in Unix timestamp
            self.lindata["Time"] = self.lindata["Time"] * 1e9  # Rectify the timestamp
        except TypeError:
            # If the data is in date format instead of Unix timestamp
            # (Old version of Neptune export)
            if self.rawlog["Timestamp"].iloc[0].day == self.rawlog["Timestamp"].iloc[-1].day:
                # If the analysis was done in 1 one day (same date)
                date = "{}-{}-{} ".format(self.rawlog["Timestamp"][0].year,
                                          self.rawlog["Timestamp"][0].month,
                                          self.rawlog["Timestamp"][0].day)

                data["year"] = date

            elif self.rawlog["Timestamp"].iloc[0].day != self.rawlog["Timestamp"].iloc[-1].day:
                # If it was done on two days
                date = "{}-{}-{} ".format(self.rawlog["Timestamp"][0].year,
                                          self.rawlog["Timestamp"][0].month,
                                          self.rawlog["Timestamp"][0].day)
                date_end = "{}-{}-{} ".format(self.rawlog["Timestamp"].iloc[-1].year,
                                              self.rawlog["Timestamp"].iloc[-1].month,
                                              self.rawlog["Timestamp"].iloc[-1].day)
                # Set the date as the first date found in the log file
                data["year"] = date
                for i in range(len(data["Time"])):
                    # Loop to find where does it go from one day to another and change the
                    # rest of the date to be like the one at the end of the log file
                    if data["Time"].iloc[i + 1][:2] == "00" and data["Time"].iloc[i][:2] == "23":
                        data.loc[i+1:, "year"] = date_end
                        break
        
            # Change the raw data format into a datetime object
            data["Time"] = pd.to_datetime(data["year"].astype(str) + " " + data["Time"].astype(str),
                                          format="%Y-%m-%d %H:%M:%S:%f")
            # Transform the dateobject into Unix timestamp
            data["Time"] = pd.to_datetime(data["Time"]).values.astype(np.int64)

            data.drop("year", axis=1, inplace=True)  # Get rid of the year column
        else:
            # If the conversion went ok, then rectify the timestamp between summer and winter time
            self.lindata["Time"] = self.lindata["Time"].astype(int)
            # Calculate the time difference between the laser log file and Neptune data
            time_difference = abs(self.lindata["Time"][0] - self.rawlog["Time"][0])
            if 3e12 < time_difference < 4e12:
                print("There is one hour difference between the laser and Neptune computer")
                self.lindata["Time"] = self.lindata["Time"] + 3.6e12 # Add one hour
            elif 6e12 < time_difference < 8e12:
                print("There is a two hour difference between the laser and Neptune computer")
                self.lindata["Time"] = self.lindata["Time"] + 7.2e12 # Add two hours
            else:
                print("No time difference between the laser and Neptune")
            
        ######################################################################
        # Now get the timelag between the laser and the Neptune and apply it.
        ######################################################################

        self.intime = (data["Time"][2]-data["Time"][1])/1e9
        print("\nThe integration time of the Neptune is: {}s\n".
              format(round(self.intime, 5)), flush=True)
        
        # Detect the standards in the analysis
        self.standard_detect()
        # Remove block ends which has wrong values
        # Every 995 integration cycles, the blocks end and a new one start
        # this operation takes between 3 and 5 seconds.
        # Removing the data point right after it is very important for lower int. times.
        td = np.append(np.diff(data["Time"]), self.intime*1e9)  # Time difference between two points (size = n-1)
        # Find where there is a >2 seconds time difference between two points,
        # as it only happens after a block ends. Then the next point is biased.
        new_ind = data["Time"].loc[td>2e9].index + 1
        # Replaced the values by NaN at the end of blocks      
        data.loc[new_ind, :] = np.ones(len(data.columns)) * np.nan 
        # Calculate the total signal had you have sharper gradient constrast
        totalsig = data.iloc[:, 1:].sum(axis=1).values
        
        # Get the location of the first Nist by calculating the minimum time
        # difference between the log entry and the Neptune time
        # (see helpers file)
        indend = helpers.getloc(data["Time"], timela_end[0])
        indstart = helpers.getloc(data["Time"], timela_start[0])
        if manstart:
            # Manual selection of the starting point of the first analysis
            # Plot the total signal of the primary standard
            ax = plt.subplot()
            x = data["Time"][indstart-50:indend+100]
            y = totalsig[indstart - 50:indend+100]
            ax.axvline(data["Time"][indstart], c="k")  # Plot where the log file
            ax.axvline(data["Time"][indend], c="k")    # recorded the start of analysis
            ax.plot(x, y)
            ax.scatter(x, y, edgecolors="k")
            my = ax.get_ylim()[1]
            ax.text(data["Time"][indstart+5], 0.9*my, "Log start" )
            ax.text(data["Time"][indend+5], 0.9*my, "Log end" )
            
            _ = Cursor(ax, useblit=True, color="k")
            coord = plt.ginput(1, timeout=-1)  # Graphical input
            plt.close()
            
            nreals = coord[0][0]
        
            timelag = data["Time"][indstart] - nreals
            
        
        else:
            # Get the location of the transition between background and signal
            # by finding the maximum gradient of the signal between start and end 
            # of the first analysis
            realstart = np.argmax(np.gradient(totalsig[indstart-10:indend+10]))
            timelag = data["Time"][realstart+indstart-10] - timela_start[0]

        # To estimate the timelag, the "realstart" has been determined ealier,
        # but it is in numpy indeces and we try to get the pandas DataFrame index
        # to match the data, so you need to apply this corrections
        
        print("The time lag between the two machines is: {}s\n".format(round(timelag/1e9, 3)),
              flush=True)

        # Correct the analysis start and end by adding the machine timelag
        anstart = timela_start - timelag
        anend = timela_end - timelag
 
        #######################################################################
        # Now separate the data into cells per analysis, with the signal,
        # background before and after.
        #######################################################################

        splitlog = {}
        self.datadi = {}
        # Loop through each analysis to create an entry in the datadi dictionary
        for i in tqdm(range(len(anstart)), desc="Sample splitting"):

            # Get the index of when the analysis start
            locstart = helpers.getloc(data["Time"], anstart[i])
            # Get the location where the sample ends
            locend = helpers.getloc(data["Time"], anend[i])

            if i == 0:
                # For the first one, include the background at the beginning

                locstartnext = helpers.getloc(data["Time"], anstart[i+1])

                self.datadi[i] = {
                    "sn" : self.sn[i],
                    "dat" : data.iloc[:locstartnext, 1:],
                    "loc" : [locstart, locend]}


            elif i == len(anstart) - 1:
                # For the last sample, just take the 100 values after the first
                # mid point, as sometimes analysis trail behind
                locendprev = helpers.getloc(data["Time"], anend[i-1])

                self.datadi[i] = {
                    "sn" : self.sn[i],
                    "dat": data.iloc[locendprev:locend+30, 1:],
                    "loc": [locstart, locend]}

            else:
                # For the rest of the data, save the slice of data and the location of
                # the sample

                locstartnext = helpers.getloc(data["Time"], anstart[i+1])
                locendprev = helpers.getloc(data["Time"], anend[i-1])
                # Save the backgrounds which extend all the way from the last
                # sample to the start of the next 
                self.datadi[i] = {
                    "sn" : self.sn[i],
                    "dat": data.iloc[locendprev+3:locstartnext-3, 1:],
                    "loc": [locstart, locend]
                    }

        # =====================================================================
        # If the check option was parsed, plot the total signal with signal, and
        # transitions highlighted
        # =====================================================================
        if check:
            # Plot the first analysis to see where is bg vs sig cuts are
            fig = plt.figure()
            time_since = data["Time"] - data["Time"][0]  # seconds since beginning of analysis
            x = time_since[self.datadi[0]["dat"].index].values # Only the first analysis
            
            plt.plot(x, self.datadi[0]["dat"]["11B"])
            # Plot the recorded signal as the indexes used later in self.bg_sub()            
            plt.axvspan(x[self.datadi[0]["loc"][0]+4], x[self.datadi[0]["loc"][1]-1],
                        color="green", alpha=0.5)
            # Same, but with the backgrounds
            plt.axvspan(x[0], x[self.datadi[0]["loc"][0]-11], color="blue", alpha=0.5)
            plt.axvspan(x[self.datadi[0]["loc"][1]+9], x[-1], color="blue", alpha=0.5)
            
            plt.xlabel("Time (s)")
            plt.ylabel("Total signal (V)")
            
            fig, ax = plt.subplots()
            ax.plot(data["Time"], totalsig)

            # Loop through each cells in the dictionary to plot the transitions
            # and the signal
            cutp = 10
            for i in self.datadi:
                start = self.datadi[i]["loc"][0]
                end = self.datadi[i]["loc"][1]

                ax.axvspan(data["Time"][start+cutp], data["Time"][end-1],
                            color="blue", alpha=0.4)
                ax.axvspan(data["Time"][start-cutp], data["Time"][start+cutp],
                            color="red", alpha=0.4)
                ax.axvspan(data["Time"][end-cutp//2], data["Time"][end+cutp//2],
                            color="red", alpha=0.4)
            
            
            plt.xlabel("Time (s)")
            plt.ylabel("Total raw voltage (v)")
            plt.title(
                "Check that all the signal, transitions and backgrounds are well defined")
            # ax2 = ax.twinx()
            # ax2.plot(data["Time"][:-1], np.diff(data["Time"]), c="k")
            plt.show()

    def bg_sub(self, mancut=False, savefig=False, endcut=1):
        """
        Background substraction.

        Take each cell defined in logtime(), remove background from the signal,
        remove outliers, and append the mean of the ratios to dataframe
        """
        # Create an empty DataFrame with all the columns to be later populated
        self.dataOvert = pd.DataFrame(columns=["sample", "N",
                                               *self.datadi[0]["dat"].columns,
                                               "11/10B", "11/10.035", "11/9.979",
                                               "11/10B2se", "11/10.0352se", "11/9.9792se"])

        ######################################################################
        #  Background substraction:
        #  Takes the signal and substract the average of surrounding backgrounds
        ######################################################################
        cutp = int(5/self.intime)  # Remove the first 5 seconds of analysis
        lens, lenb = [], []
        sig_bg_ratio = np.zeros(len(self.datadi))

        for i in tqdm(self.datadi, desc="Background substraction"):

            if mancut:
                """ 
                Manual cutting of the signal by displaying samples one by one
                with the raw 11/10B ratio and 11B signal                    
                """
                # Extract the signal from the automatic cut
                sig = self.datadi[i]["dat"].loc[self.datadi[i]["loc"][0]+cutp:
                                                self.datadi[i]["loc"][1]-endcut]

                # Plot the raw ratios
                rr = sig["11B"] / sig["10B"]  # Calculate the raw ratio
                    
                # TODO Downhole fractionation of low [B] samples?
                
                sn = self.datadi[i]["sn"]
                # Plot the single analysis
                fig = helpers.an_plot(sig, i, len(self.datadi), sn)
                
                plt.tight_layout()  # Magic to look nice
                # Use an interactive cursor to select the start and end
                _ = Cursor(fig.axes[-1], useblit=True, color="k", linewidth=1)
                coord = plt.ginput(2, timeout=-1)  # Get the coords
                
                if coord: # If points were selected    
                    if len(coord) == 1:
                        start, end = 0, round(coord[0][0])
                    else:
                        start, end = round(coord[0][0]), round(coord[1][0])
                        if start < 0: # If misclicked below 0 rectify.
                            start = 0
                            
                    sig = sig.iloc[start:end+1, :]
                    
                if savefig:
                    try:
                        ndir = self.directory + "/man_cut_plots/"
                        os.mkdir(ndir)
                    except FileExistsError:
                        pass
                    
                    if coord:    
                        plt.axvline(start, ls="--", c="k")
                        plt.axvline(end, ls="--", c="k")
                    else:            
                        plt.axvline(0, ls="--", c="k")
                        plt.axvline(len(rr), ls="--", c="k")
                
                
                    plt.savefig(ndir + f"{sn}_{i}.png")
                    
                plt.close(fig)

            else: # Automatic cutting
                # Extract the signal from the whole sample and remove the edges
                sig = self.datadi[i]["dat"].loc[self.datadi[i]["loc"][0]+cutp:
                                                self.datadi[i]["loc"][1]-endcut]
                # If the reported value equals exactly 0, which is normally impossible
                # since there is always electronic noise, it means that the detector
                # has been saturated, skip this analysis. 
                sat = sig["11B"].loc[sig["11B"] == 0.0]
                if sat.count():
                    print("The faradays were saturated, skipped")
                    empt_row = [self.sn[i]] + (self.dataOvert.shape[1]-1) * [np.nan]
                     
                    self.dataOvert.loc[len(self.dataOvert)] = empt_row
                    continue
                
                if savefig:
                    try: # Save the figures in new directory
                        ndir = self.directory + "/auto_cut_plots/"
                        os.mkdir(ndir)
                    except FileExistsError:
                        pass # If it already exists, skip
                    
                    sn = self.datadi[i]["sn"]
                    # Plot the single analysis
                    fig = helpers.an_plot(sig, i, len(self.datadi), sn)
                    
                    plt.tight_layout()  # Magic to look nice
                    plt.savefig(ndir + f"{sn}_{i}.png")
                    plt.close()
                    
            # !!! Changed the tail n rows to be only half of the front cut
            lens.append(len(sig))
            # Extract the backgrounds before and after the signal
            bgs = [self.datadi[i]["dat"].loc[:self.datadi[i]["loc"][0]-11],
                   self.datadi[i]["dat"].loc[self.datadi[i]["loc"][1]+9:]]
            # Trim the backgrunds to match David version
            bgs = [bgs[0].iloc[6:, :], bgs[1].iloc[:-8, :]]
            
            if len(sig) == 0:
                """
                If the cell is empty, there must be problem with the log file
                extending beyond the data of the Neptune, just raise warning
                """
                warnings.warn(f"Sample {self.datadi[i]['sn']} appears to have 0"
                              " data points, perhaps an error on the log file?")
                self.datadi[i]["sig"] = []
                self.dataOvert.loc[len(self.dataOvert)] = [self.datadi[i]["sn"]]+[0]*(len(self.dataOvert.columns) - 1)

            else:
                #! Difference between the two backgrounds is here
                # The problem is that in the matlab script, when there is one
                # outlier, the whole row is not counted (all cups).
                # bg_removal = [helpers.matsdrem(bg, 2) for bg in bgs]
                # Remove outliers in background. Calculated per cup
                bg_removal = [helpers.stdmask(bg, 2) for bg in bgs]  
                # Average of the backgrounds 
                avgbg = [x.mean() for x in bg_removal]
                # Loop through the channels and average the two backgrounds
                bgval = [np.mean([avgbg[0][i], avgbg[1][i]])
                         for i in range(len(avgbg[0]))]
                # This should be the same treatment as the background.
                # sig = helpers.stdmask(sig, 3)   # Remove 3 sd in signal
                sig = helpers.matsdrem(sig, 3)  # David's sd removal
                
                sig -= bgval   # Background substraction

                # Make the ratios based on the analyte
                try:
                    ratiobca = sig["11B"]/sig["10.035"]
                    ratiobca9 = sig["11B"]/sig["9.979"]
                except KeyError:  # If it is an old file with old cup config
                    ratiobca = sig["11B"]/sig["10.084"]
                    ratiobca9 = sig["11B"]/sig["9.95"]
                    
                ratio = sig["11B"] / sig["10B"]
                # ratio_fil = helpers.stderem(ratio, 2)
                # This second is the 2sd outlier removal, the first is the
                # removal of spikes
                # ratio.loc[abs(ratio - ratio.mean()) > ratio.std()*2] = np.nan

                b2se = ratio.sem()*2  # 2SE of the 11/10 ratio
                # 2 stantdard error of the 11B/10.035 and 11B/9.979 ratio
                bca2se = ratiobca.sem()*2
                bca92se = ratiobca9.sem()*2
                # Add the filtered ratio values to the "sig" Series
                sig.insert(len(sig.columns), "11/10B", ratio)
                sig.insert(len(sig.columns), "11/10.035", ratiobca)
                sig.insert(len(sig.columns), "11/9.979", ratiobca9)
            
                new_row = [self.datadi[i]["sn"], len(sig), *sig.mean(), b2se, bca2se, bca92se]
                # Append the sample name, means of all channels and 2se to dataframe
                self.dataOvert.loc[len(self.dataOvert)] = new_row
                # Calculate 11B signal to background ratio
                sig_bg_ratio[i] = sig["11B"].mean() / bgval[-1]
                
                self.datadi[i]["sig"] = sig
                self.datadi[i]["bg"] = bgval
                self.datadi[i]["len"] = len(sig)
                
        # Insert the laser parametres
        self.dataOvert.insert(1, "spotsize", self.spotsizes)
        self.dataOvert.insert(2, "fluence", self.fluences)
        self.dataOvert.insert(3, "rep_rate", self.rep_rates)
        self.dataOvert.insert(4, "sig_bg", sig_bg_ratio)
        ###
        # Plot backgrounds
        ###
        labels = ["9.979", r"$^{10}$B","10.035","10.433",r"$^{11}$B"]
        bg_vals = np.array([self.datadi[i]["bg"] for i in self.datadi])

        fig, ax = plt.subplots()
        for i in range(len(labels)):
            ax.plot(bg_vals[:, i]*1e3, label=labels[i])
            
        ax.set_ylabel(r"Background (mV)")
        ax.set_xlabel("Analyte #")
        ax.set_xlim(0, len(self.datadi))
        # ax.set_yscale("log")
        ax.legend()
        plt.savefig(self.directory + "/background_drift.png")

    def nistcor(self, save=True, poly=False):
        """
        Correct for drift and mass bias using the NISTs brackets.

        For each spot size, locate where the primary standard is (usually NIST612)
        and create nist blocks. The blocks are then averaged and used to correct
        the samples between the two brackets.
        """
        # TODO Record fractination factor. As f of time ? 
        # TODO calculate sensitivity of NIST, as f of time.
        nistloc = self.dataOvert[self.dataOvert["sample"].str.match("nist")].index
        # Transfer the 11/10 ratios into a numpy array because it is easier to deal with
        bvals = self.dataOvert["11/10B"].values
        b11sig = self.dataOvert["11B"].values  # Do the same with the raw 11B signal
        corrected = np.zeros(len(bvals))  # Empty NIST corrected values
        bconcs = np.zeros(len(bvals))     # Empty array for [B]
        c_err = np.zeros(len(bvals))
        NISTratio = 4.041648165  # The accepted 11/10 of NIST612
        
        # =====================================================================
        # Outlier removal of the nists        
        # =====================================================================
        # TODO NIST trend line? 

        nistvals = self.dataOvert["11/10B"][nistloc].values
           
        outlier = True
        
        if outlier:
            zs = abs(scipy.stats.zscore(nistvals))  # Calculate the zscore of the nists
            if (zs > 3).any():  # If there are 2SD outliers
                print(f"There are {len(zs[zs > 2])} NIST outliers")
                outloc = np.where(zs>2)[0]
                nistloc = np.delete(nistloc, outloc)
                nistvals = np.delete(nistvals, outloc)
            else:
                print("No NIST outliers")
        

        # =====================================================================
        # Group the NISTs brackets together
        # =====================================================================

        nistblocks = []
        block = []

        for i in range(len(nistloc)):
            if i == 0:  # For the first one only
                block.append(nistloc[i])
                continue

            if nistloc[i] - nistloc[i-1] == 1:  # If the current nist has a neighbour
                block.append(nistloc[i])
                if i == len(nistloc) - 1:  # If it is the last one
                    nistblocks.append(block)
            else:
                # When the difference of location is more than 1
                nistblocks.append(block)
                block = [nistloc[i]]      # i.e. the next NIST bracket

                if i == len(nistloc)-1:  # If the last one is alone, just append it
                    nistblocks.append(block)

        # =====================================================================
        # NIST bracket correction
        # =====================================================================
        nist_fr = []
        for i in range(len(nistblocks)-1):

            nb = np.mean(bvals[nistblocks[i]])  # NISTs before
            na = np.mean(bvals[nistblocks[i+1]])  # NISTs after
            nm = np.mean((na, nb))  # NIST means
            frac_factor = NISTratio/nm
            nist_fr.append(frac_factor)
            # Take all the values between the two NISTs brackets and apply the
            # mass bias correction to it
            newvals = bvals[nistblocks[i][0]:nistblocks[i+1][-1]+1] * frac_factor
            # Fill the empty "corrected" array with the correct values
            corrected[nistblocks[i][0]:nistblocks[i+1][-1]+1] = newvals
            
            # Here, make an average of the two NIST pairs 11B signal
            nm11 = np.mean((np.mean(b11sig[nistblocks[i]]),
                           np.mean(b11sig[nistblocks[i+1]])))    
            
            # Make the ratio of the 11B signal of the sample to the NIST average
            bcest = (b11sig[nistblocks[i][-1]+1:nistblocks[i+1][0]]/nm11)
                       
            # Fill the array with the above calculated ratio times the 
            # boron concentration in NIST612
            bconcs[nistblocks[i][-1]+1:nistblocks[i+1][0]] = bcest * 34.3


        # Put the NIST values back.
        # Since the correction only applies to what is in between the brackets,
        # The "corrected" array has 0 where the NISTS are.
        # corrected[np.where(corrected == 0)] = bvals[np.where(corrected == 0)]
        # =============================================================================
        # Polynomial fitting through the nist trial -> like Jie        
        # =============================================================================
        order = 5  # Order of the polynomial
        popt, pcov = scipy.optimize.curve_fit(helpers.polynomfit, nistloc, bvals[nistloc],
                                              p0=[0]*order)
        
        nist_pol = helpers.polynomfit(np.array(self.dataOvert.index),*popt)
        nist_cor_fac = NISTratio/nist_pol
        # TODO the same with [B]
        if poly:
            self.dataOvert["11/10Bc"] = bvals * nist_cor_fac
        else:
            self.dataOvert["11/10Bc"] = corrected
            
        self.dataOvert["[B]"] = bconcs
        # Transform into delta notation using NISTRM 951
        self.dataOvert["d11B"] = (self.dataOvert["11/10Bc"]/4.04367 - 1)*1000
        # Our way of calculating SE in permil is 1:1 with Jie's
        self.dataOvert["r2se"] = (self.dataOvert["11/10B2se"]/self.dataOvert["11/10B"])*1000
        # =============================================================================
        # B/Ca calculate by carbonate polynomial fitting through session
        # =============================================================================
        
        std_bca_values = {"jct": unc.ufloat(191, 9),  # umol/mol from Hathorne 2013 G3, 2SD
                          "jcp": unc.ufloat(460, 23),  # umol/mol from Hathorne 2013 G3, 2SD
                          "macs": unc.ufloat(56, 2)}  # umol/mol from Standish 2019 RCMS, 2SD
        # Create a separate DataFrame for B/Ca, too many columns, will render the main
        # df too crowded.
        bca_res_df = pd.DataFrame()
        bca_res_df["sample"] = self.dataOvert["sample"]
        if poly: # Use a polynomial fit through the standards 10.035
           for i, std in enumerate(std_bca_values): # Loop through each standard
                std_val = self.dataOvert[self.dataOvert["sample"].str.match(std)]
                std_loc = std_val.index
                order = 5  # Order of the polynomial
                popt, pcov = scipy.optimize.curve_fit(helpers.polynomfit, std_loc, std_val["11/10.035"],
                                                      p0=[0]*order)
                session_pol = helpers.polynomfit(np.array(self.dataOvert.index),*popt)
                bca_fac = std_bca_values[std] / session_pol
                b11_10u = unp.uarray(self.dataOvert["11/10.035"], self.dataOvert["11/10.0352se"])
                bca_pol = b11_10u * bca_fac
                bca_poln = unp.nominal_values(bca_pol)
                bca_polu = unp.std_devs(bca_pol)
                bca_res_df["B/Ca_pol_" + std] = bca_poln
                bca_res_df["B/Ca_pol_" + std + "_2se"] = bca_polu
        
        # =============================================================================
        # B/Ca calculation by carbonate standard bracketing
        # =============================================================================
        # B/Ca estimates
        # Follow Babila et al. 2021 (Sci. Adv.) by bracketing JCt
        for i, std in enumerate(std_bca_values):
            std_loc = self.dataOvert[self.dataOvert["sample"].str.match(std)].index
            casig = self.dataOvert["11/10.035"].values   # Extract 11/10.035 values
            casigse = self.dataOvert["11/10.0352se"].values   # Extract 11/10.035 values 2SE
            bcas = np.zeros_like(bvals)   # Empty B/Ca array
            bcau = np.zeros_like(bvals)   # Empty B/Ca uncertainty array
        
            for i in range(len(std_loc) - 1):
                # B/Ca values in between the two JCt
                # bcasamp = casig[std_loc[i]+1:std_loc[i+1]]
                # bcasesamp = casigse[std_loc[i]+1:std_loc[i+1]]
                bca_u = unp.uarray(casig[std_loc[i]+1:std_loc[i+1]],        # 11/10.035 nominal value  
                                   casigse[std_loc[i]+1:std_loc[i+1]])      # 11/10.035 std devs
                # B/Ca values of the first and second JCt of the bracket
                bca_std_1 = casig[std_loc[i]]
                bca_std_2 = casig[std_loc[i+1]]

                bca_jct_val = np.average((bca_std_1, bca_std_2))
                    
                stdval = std_bca_values[std]
                bca_values = (bca_u/ bca_jct_val) * stdval
                bcas[std_loc[i]+1:std_loc[i+1]] = unp.nominal_values(bca_values)
                bcau[std_loc[i]+1:std_loc[i+1]] = unp.std_devs(bca_values)
            
            bca_res_df["B/Ca_bra_" + std] = bcas
            bca_res_df["B/Ca_bra_" + std + "_2se"] = bcau
        
        self.dataOvert["B/Ca"] = bcas # Save results to the main dataset.
        self.dataOvert["B/Ca_2sd"] = bcau

        # 2nd method for calculating B/Ca, use the [B] calibrated using NIST612
        # Because NIST612 is more homogeneous in terms of []
        # Calculate offset between [B] measured and reported [B] in JCt
        # Apply this offset either in bracketed form or average all and use it as is
        # jct [B] = 23
        # jcp [B] = 50 ug/g
        # macs3 [B] = 8.2
        # jct = self.dataOvert.loc[self.dataOvert["sample"].str.match("jct")]
        # jctoff = jct["[B]"] - 23 
        
        # =============================================================================
        #         # Plot the changes of nist values through the analysis.
        # =============================================================================
        x = nistloc
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        # plt.scatter(x, nistvals, ec="k")
        ax[0].errorbar(x, corrected[nistloc], self.dataOvert.loc[nistloc, "r2se"].values/1e3,
                       mec="k", fmt="o")
        
        # Plot a regression line
        m, c = np.polyfit(x, corrected[nistloc], 1)
        ax[0].plot(x, x*m + c, c="k")
        # Plot the mean and 2sd of the whole nists analysis
        nm = corrected[nistloc].mean()
        nsd = corrected[nistloc].std()*2
        ax[0].axhline(nm, ls="--", color="black")
        ax[0].axhline(nm + nsd, ls=":", color="black")
        ax[0].axhline(nm - nsd, ls=":", color="black")
        
        ax[0].set_title(r"Corrected NIST612 $^{11}$B/$^{10}$B variability")
        # plt.ylabel(r"relative $\delta ^{11}$B of NIST612")
        ax[0].set_ylabel("Corrected 11/10 B of NIST612")
        ax[0].set_xlabel("Analysis #")
        labs = [item.get_text() for item in ax[0].get_yticklabels()]
        ax[0].set_ylim(ax[0].get_ylim())
        ax2 = ax[0].twinx()
        
        ax2.set_ylim(ax[0].get_ylim())
        ax2.set_ylabel(r"$\delta^{11}$B$_{951}$")
        
        nlabs = (np.array(ax[0].get_yticks())[1:-1]/4.04367 - 1)*1000
        nlabs = [f"{lab:0.1f}" for lab in nlabs]
        ax2.set_yticks(ax[0].get_yticks()[1:-1])
        ax2.set_yticklabels(nlabs)
        ###########
        # Plot the signal intensity
        ###########
        ax[1].plot(x, b11sig[nistloc]*1e3, "o", mec="k")
        # ax[1].plot(x, (x*ols_nist_11B.params[1] + ols_nist_11B.params[0])*1e3, c="k")
        # ax[1].plot(px, pfun(px, *B11popt)*1e3, "r")
        
        ax[1].set_xlabel("Analysis #")
        ax[1].set_ylabel(r"NIST612 $^{11}$B signal (mV)")
        ax[1].set_title(r"NIST612 $^{11}$B signal through session")
        fig.tight_layout()
        if save:
            plt.savefig(self.directory + "/nist_drift_plot.png")
            
        
    def conc_plot(self, save=False):
        """Plot B concentration of all the samples."""
        plt.figure()
        # plt.scatter(self.dataOvert["11/10.035"], self.dataOvert["[B]"], ec="k",
                    # label="Samples", alpha=0.7)
        plt.scatter(self.dataOvert["11/10.035"], self.dataOvert["11B"], ec="k",
                    label="Samples", alpha=0.7)
                
        for s in self.stds:
            if s == "nist612":
                continue
            std_data = self.dataOvert.loc[self.dataOvert["sample"].str.match(s)]
            # plt.scatter(std_data["11/10.035"], std_data["[B]"], ec="k", label=s, s=45)
            plt.scatter(std_data["11/10.035"], std_data["11B"], ec="k", label=s, s=45)
            
        plt.legend()
        plt.xlabel(r"11/10.035 $\propto$ B/Ca")
        plt.ylabel(r"$^{11}$B (v)")
        plt.xlim(0,1800)
        plt.savefig("Boron_concentration_plot.pdf")
        plt.show()
        

    def ca_cor(self, outlier=True, save=True, manout=False, ca9=False, yel=False,
               spot=False, flu=False):
        """
        Correction of the Ca matrix effect.

        Correct for the Matrix effect using the methodology from
        Evans et al. 2021
        """        
        # Make the 'B/Ca' ratio or 'signal/Ca interference'
        # bca = self.dataOvert["11B"]/self.dataOvert[ca]
        # bca2 = self.dataOvert["11B"]/self.dataOvert[ca_ch[0]]
        # bca = self.dataOvert["11B"]
        bca = self.dataOvert["11/10.035"]
        bca2se = self.dataOvert["11/10.0352se"]
        bca9 = self.dataOvert["11/9.979"]
        bca92se = self.dataOvert["11/9.9792se"]
        # bca9 = self.dataOvert["11B"] / self.dataOvert["9.95"]
        # Uncertainty in bca
        # Dictionary of the solution values for the standards
        stds = {"jcp": 24.3,
                "jcpnp": 24.3,
                "jct": 16.3,
                "jctnp": 16.3,
                "macs-3": 0.57,
                "macs3": 0.57,
                "mac3np": 0.57,
                "macs-3 nano": 0.57,
                "blue": -0.12}
                # "yellow":-18.7}
        
        stds = {"jcp": [24.36, 0.055],  # Mean +- 1SD
                "jct": [16.3, 0.30],    # Mean +- 1SD
                "macs": [0.57, 0.055],
                }
        if yel:
            stds["yellow"] = [-18.7, 0.2]
            # "blue": [-0.12, 0.41]
        
        mx, my = [], []
        mx2se, my2se = [], []
        mx9, mx92se = [], []
        stdfound = []
        std_locs = {}
        
        if spot:
            # Filter the dataset for a specific spot size
            dataset = self.dataOvert[self.dataOvert["spotsize"] == spot]
        elif flu:
            # Filter the dataset for a specific fluence
            dataset = self.dataOvert[self.dataOvert["fluence"].between(flu-1, flu+1)]
        else:
            dataset = self.dataOvert
        # Loop through the standards to find their location in the sample list
        # and save their B/Ca and d11B values to later fit the model
        for i in stds:
            # loc = np.where(dataset["sample"].str.match(i))[0]
            loc = dataset[dataset["sample"].str.match(i)].index
            
            if len(loc) == 0:
                print(f"Did not find the standard: {stds[i]}")
                continue  # If did not find the standard ignore
            
            print(f"{i} standard:")
            # Remove any sample that has negative 11/10.035 (or 9.979)
            above_zero = loc[np.where(bca[loc].values > 0)]
            print("{} samples with 11/10.035 < 0".format(len(loc) - len(above_zero)))
            # Remove any sample that has negative SE of d11B
            pos_se = above_zero[np.where(dataset.loc[above_zero, "r2se"] > 0)]
            print("{} samples with negative uncertainty".format(len(above_zero) - len(pos_se)))
            # Remove any sample with a Dd11B > 0 
            out_loc = pos_se[np.where(dataset.loc[pos_se, "d11B"] - stds[i][0] < 0)]
            print("{} samples with Dd11B > 0\n".format(len(pos_se) - len(out_loc)))
            std_locs[i] = out_loc
            
            mx.append(bca[out_loc].values)
            mx2se.append(bca2se[out_loc].values)
            mx9.append(bca9[out_loc].values)
            mx92se.append(bca92se[out_loc].values)
            
            uy1 = unp.uarray(dataset["d11B"][out_loc], dataset["r2se"][out_loc]/2)
            uy2 = unp.uarray(stds[i][0], stds[i][1])
            uy = uy1 - uy2
            myfil = unp.nominal_values(uy)
            my.append(myfil)
            my2sefil = unp.std_devs(uy)*2
            my2se.append(my2sefil)
            
            stdfound.append(i)
            
        # Concatenate all three standards into a single array.
        mxc = np.concatenate(mx)
        mxc2se = np.concatenate(mx2se)
        mxc9 = np.concatenate(mx9)
        mxc92se = np.concatenate(mx92se)
        
        myc = np.concatenate(my)
        myc2se = np.concatenate(my2se)
        # Use the 9.979 channel for correction instead of 10.035
        if ca9:
            mxc = mxc9
            mx = mx9
            mx2se = mx92se
            
        if len(mx) == 0 or len(my) == 0:
            raise ValueError("There are no Jcp, Jct or MACS-3 in this "
                             "spot size, so there is no correction.")

        if manout:
            print("Manual outlier removal")
            fig, ax = plt.subplots()
            
            for i in range(len(mx)):
                ax.scatter(mx[i], my[i], label=stdfound[i])
                # ax.scatter(mx[i], my[i], ec="k", s=ss[i]*2, alpha=0.6, label=stdfound[i])
                # ax.errorbar(mx[i], my[i], my2se[i], mx2se[i], mec="k",
                            # fmt="o", label=stdfound[i])
            
            if ca9:
                plt.xlabel(r"B/'Ca' ($^{11}$B/9.979)")
            else:
                plt.xlabel(r"B/'Ca' ($^{11}$B/10.035)")
            plt.ylabel(r"$\delta ^{11}$B")
            plt.legend()
            # Plot the data and interactively choose points
            _ = Cursor(ax, useblit=True, color="k")
            coord = plt.ginput(-1, timeout=-1)  # Graphical input
            plt.close(fig)
            # Stack the x and y values to pass it to a cKDTree
            xy = np.column_stack([mxc, myc])
            xy2se = np.column_stack([mxc2se, myc2se])
            ind = []
            out = []
            # If points were selected, then loop over each and delete them
            if coord:
                for j in coord:
                    dists= np.sum( (xy-j)**2, axis=1)
                    minind = np.argmin(dists)
                    ind.append(minind)
                    out.append(xy[minind]) # Save xy values of outliers
            
            out = np.array(out)
        
            xy = np.delete(xy, ind, 0)  # Delete the outliers
            xy2se = np.delete(xy2se, ind, 0) # Delete outliers in 2se array
            # Destack the arrays
            mxc = xy[:, 0]
            myc = xy[:, 1]
            mxc2se = xy2se[:, 0]
            myc2se = xy2se[:, 1]

        try:
            ### Fit the non-linear regression with the filtered data
            popt, pcov = scipy.optimize.curve_fit(helpers.exp_fun, mxc, myc, p0=[-2000, -1, 1],
                                                  absolute_sigma=True, maxfev=int(1e6))
        except RuntimeError:
            warnings.warn(f"Could not converge the non linear regression."\
                  " This might be due to a great degree of scatter / outliers,"\
                  " try using the outlier removal.")
            lens = len(self.dataOvert["11B"])
            self.dataOvert[f"11/10.035"] = np.ones(lens)*np.nan
            self.dataOvert["Dd11B"] = np.ones(lens)*np.nan
            self.dataOvert["cald11B"] = np.ones(lens)*np.nan
            self.dataOvert["r2se"] = np.ones(lens)*np.nan
        

        multiple_nl_reg = False
        if multiple_nl_reg:
            # Try to consider both the 9.979 and 10.035 cup in consideration for 
            # the regression. It is not working too well, large d11B difference
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.ensemble import GradientBoostingRegressor
            from sklearn import tree

            uwcloc = self.dataOvert[self.dataOvert["sample"].str.match("uwc-3")].index 
            # blueloc = self.dataOvert[self.dataOvert["sample"].str.match("blue")].index 
            # sec_std_loc = np.concatenate((uwcloc, blueloc))
            xc = self.dataOvert.loc[uwcloc,"11/10.035"]
            x_test = self.dataOvert.loc[uwcloc,["11/9.979","11/10.035"]]
            y_test = self.dataOvert.loc[uwcloc, "d11B"] - 20.25
            
            xd = pd.DataFrame(list(zip(mxc9, mxc)), columns=["11/9.979", "11/10.035"])
            plt.figure()
            plt.plot(xd.iloc[:, 0] - xd.iloc[:, 1], "o")
            clf = RandomForestRegressor(max_depth=10, n_estimators=1000)
            clf.fit(xd, myc)
            
            # clfres = myc - clf.predict(xd)
            clfres = y_test - clf.predict(x_test)
            print(clf.score(xd, myc))   
            
            gbr = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1,
                                            max_depth=5, random_state=0,loss='squared_error')
            gbr.fit(xd, myc)
            # gbrres = myc - gbr.predict(xd)
            gbrres = y_test - gbr.predict(x_test)
            print(gbrres)
            print(gbr.score(xd, myc))
            
            dtr = tree.DecisionTreeRegressor()
            dtr.fit(xd, myc)
            # dtrres = myc - dtr.predict(xd)
            dtrres = y_test - dtr.predict(x_test)
            print(dtr.score(xd, myc))

            pred = helpers.exp_fun(mxc, *popt)
            res = myc - pred
            plt.figure()
            
            # plt.plot(mxc, res, "o", mec="k", label="curve fit")
            plt.plot(xc, clfres, "o", mec="k", label="Random Forest")
            plt.plot(xc, gbrres, "o", mec="k", label="Gradient Boost")
            plt.plot(xc, dtrres, "o", mec="k", label="Decision Tree")
            plt.axhline(0, ls="--", c="k")
            plt.legend()
            
            ax = plt.figure().add_subplot(projection='3d')
            ax.scatter(xd.iloc[:, 0], xd.iloc[:, 1], myc)
            
        mc_cal = False
        if mc_cal:
            # Monte Carlo approach to incorporate the uncertainty of Dd11B and B/Ca
            n = int(1e4)
            mc_xr = []
            mc_yr = []
            for i in range(len(mxc)):
                mc_xr.append(np.random.normal(mxc[i], mxc2se[i]/2, n))
                # mc_xr.append(np.random.normal(mxc[i], 0, n))
                mc_yr.append(np.random.normal(myc[i], myc2se[i]/2 , n))
                # mc_yr.append(np.random.normal(myc[i], 0, n))
            mc_xr = np.array(mc_xr)
            mc_yr = np.array(mc_yr)

            coefs = np.zeros((n, 3))
            covs = np.zeros((n,3,3))
            # m, c = np.zeros(n), np.zeros(n)  
            px = np.linspace(min(mxc), max(mxc), 100)
    
            for i in tqdm(range(n)):
                coefs[i, :], covs[i, :, :] = scipy.optimize.curve_fit(helpers.exp_fun, 
                                                                      mc_xr[:, i],
                                                                      mc_yr[:, i],
                                                                      p0=popt,
                                                                      absolute_sigma=True,
                                                                      maxfev=10000)
            
            mc_popt = np.percentile(coefs, 50, axis=0)
            mc_pcov = np.percentile(covs, 50, axis=0)
            
            
            
            mc_pars = unc.correlated_values(mc_popt, mc_pcov)
            mc_n_pars = unp.nominal_values(mc_pars)
            # Get the parametres with their uncertainties
            mcnpars = unp.nominal_values(mc_pars)
            mcpy = helpers.exp_fun(px, *mc_pars)
            nommc = unp.nominal_values(mcpy)
            stdmc = unp.std_devs(mcpy)
        
        ########################################
        pars = unc.correlated_values(popt, pcov)
        npars = unp.nominal_values(pars)
        
        # Generate some x values to plot the regression curve
        px = np.linspace(min(mxc), max(mxc), 1000)
        # =============================================================================
        # Fit the artificial x values to the model with uncertainties in the parameters
        # =============================================================================
        py = helpers.exp_fun(px, *pars)
        nom = unp.nominal_values(py)
        std = unp.std_devs(py)    
        
        if outlier:
            # Automatically remove outliers
            # Calculate residuals
            pred = helpers.exp_fun(mxc, *npars)
            res = myc - pred
            mres = []
            # Values to be fitted: Average of 11B/10.035 vs average 2SE for that range of B/Ca values
            fx, fy = np.zeros(len(stds)), np.zeros(len(stds))
            for i, v in enumerate(stds):
                fx[i] = np.mean(mx[i])
                loc = std_locs[v]
                fy[i] = self.dataOvert.loc[loc, "r2se"].mean()
            
            # Fit the an exponential curve through the average uncertainty for a specific
            # standard. The lower the B/Ca, the higher the d11B uncertainty.
            res_popt, res_pcov = scipy.optimize.curve_fit(helpers.exp_fun, fx, fy, p0=[-2000, -1, -1],
                                                          absolute_sigma=True, maxfev=int(1e6))
            # The fit was 
            out_v = helpers.exp_fun(mxc, *res_popt)
            # Filter out the 2SD of the fit
            fil_mask = abs(res)<out_v*2
            out_mask = abs(res)>out_v*2
            
            mxc_fil = mxc[fil_mask]
            res_fil = res[fil_mask]
            mxc_out = mxc[out_mask]
            res_out = res[out_mask]
            
            # Plot the residuals and the outliers
            plt.figure()
            res_px = np.linspace(min(mxc), max(mxc))
            for i in range(len(mx)):
                pred = helpers.exp_fun(mx[i], *npars)
                res = my[i] - pred
                mres.append(np.mean(res))
                plt.errorbar(mx[i], res, my2se[i], mx2se[i], fmt="none", alpha=0.6, c="gray", zorder=1)
                plt.scatter(mx[i], res, ec="k", zorder=2)
                plt.axhline(0, ls="--", c="k")
            
            # Plot the outliers on top
            plt.plot(mxc_out, res_out, "ro")
            # plot the residual SD and 2SD
            plt.plot(res_px, helpers.exp_fun(res_px, *res_popt), c="r", zorder=3, ls="--")
            plt.plot(res_px, -helpers.exp_fun(res_px, *res_popt), c="r", zorder=3, ls="--")
            plt.plot(res_px, 2*helpers.exp_fun(res_px, *res_popt), c="r", zorder=3, ls=":")
            plt.plot(res_px, -2*helpers.exp_fun(res_px, *res_popt), c="r", zorder=3, ls=":")
            # Axis labelling
            plt.ylabel(r"$\Delta \delta ^{11}$B residuals")
            plt.xlabel(r"$^{11}$B/10.035")
            plt.title("Outlier removal based on residuals from the fit")
            if save:
                plt.savefig(f"{self.directory}/Ca_cor_outlier_removal.png", dpi=200)
            ########
            # Refit the filtered data
            mxc = mxc[fil_mask]
            myc = myc[fil_mask]
        
        popt, pcov = scipy.optimize.curve_fit(helpers.exp_fun, mxc, myc, p0=[-2000, -1, 1],
                                              absolute_sigma=True, maxfev=10000)
        pars = unc.correlated_values(popt, pcov)
        npars = unp.nominal_values(pars)
        upars = unp.std_devs(pars)
        
        # plot data 
        fig, ax = plt.subplots()
        for i in range(len(mx)):
            plt.errorbar(mx[i], my[i], my2se[i], mx2se[i], fmt="none", color="k", alpha=0.6)
            plt.scatter(mx[i], my[i], label=stdfound[i], zorder=2, edgecolor="k")
        
        if manout:
            if out.any():
                plt.scatter(out[:, 0], out[:, 1], facecolor="r", ec="r", s=50, zorder=3)
        # Calculate the prediction bands 95% confidence
        lpb, upb = helpers.predband(px, mxc, myc, npars, helpers.exp_fun)

        # plot the regression
        plt.plot(px, nom, c='black', label="normal")
        if mc_cal:
            plt.plot(px, nommc, label="50th percentile")
        # plt.plot(px, nommcm, label="mean")
        plt.legend()
        plt.ylabel(r"$\Delta \delta ^{11}$B")
        plt.xlabel(r"$^{11}$B/10.084 or ~'B/Ca'")        
        # uncertainty lines (95% confidence)
        plt.plot(px, nom - 2*std, c='orange',
                 label='95% Confidence Region')
        plt.plot(px, nom + 2* std, c='orange')
        # prediction band (95% confidence)
        plt.plot(px, lpb, 'k--', label='95% Prediction Band', alpha=0.5)
        plt.plot(px, upb, 'k--', alpha=0.5)
        
        # Plotting details
        plt.ylabel(r"$\Delta \delta ^{11}$B")
        plt.xlabel(r"$^{11}$B/10.084 or ~'B/Ca'")
        plt.legend(loc='best')
        plt.grid(False)

        # Save the difference of d11B compared to true by fitting the raw
        # d11B values to the non-linear regression determined above
        if mc_cal:
            self.dataOvert["Dd11B"] = unp.nominal_values(
                helpers.exp_fun(self.dataOvert["11/10.035"], *mc_n_pars))
            
        else:
            if ca9: # Use the 9.979 mass
                self.dataOvert["Dd11B"] = unp.nominal_values(
                    helpers.exp_fun(self.dataOvert["11/9.979"], *npars))
            else: # Use the 10.035 mass
                self.dataOvert["Dd11B"] = unp.nominal_values(
                    helpers.exp_fun(self.dataOvert["11/10.035"], *npars))
    
        if multiple_nl_reg:
            Dd11bm = clf.predict(self.dataOvert[["11/9.979", "11/10.035"]])
            Dd11bm = gbr.predict(self.dataOvert[["11/9.979", "11/10.035"]])
            self.dataOvert["Dd11Bm"] = Dd11bm
            self.dataOvert["cald11Bm"] = self.dataOvert["d11B"] - Dd11bm
            
            
        # Save the calibration d11B by substracting the difference from true
        self.dataOvert["cald11B"] = self.dataOvert["d11B"] - self.dataOvert["Dd11B"]
        
        # Save the regression parameters
        reg_pop = pd.DataFrame(list(zip(npars, upars)), columns=["nom", "sd"],
                               index=["a", "b", "c"])
        reg_pop.to_csv(f"{self.directory}/Dd11B_regression_coefs.csv")
        
        if save:
            if yel:
                plt.savefig(f"{self.directory}/Ca_cor_yellow.png")
            else:
                plt.savefig(f"{self.directory}/Ca_cor.png", dpi=200)
            
    def samp_plot(self, name, save=False):
        """Display a single analysis, by name."""
        cou = self.sn.count(name)
        if cou:
            if cou == 1:
                ind = self.sn.index(name)
                dat = self.datadi[ind]
                sig = dat["sig"]
                sn = dat["sn"]
                # Plot the single analysis
                fig = helpers.an_plot(sig, ind, len(self.datadi), sn)
                
                plt.tight_layout()  # Magic to look nice
            else:
                print("More than one samples has this name")
                # TODO Plot all the samples with this name
            
        else:
            print("Sample not in list")

    def uwc_plot(self, save=True, name=False):
        """Plot UWC-1 and 3, and their solution values."""
        # Extract the uwcs values and calculate their mean and 2se
        # UWC-3
        uwc3strings = ["uwc-3", "uwc3", "uwc_3"]
        for string in uwc3strings:  # Find the samples in spite of the mispellings
            uwc3 = self.dataOvert.loc[self.dataOvert["sample"].str.match(string)]
            if not uwc3.empty:
                break
            
        uwc3vals = uwc3["cald11B"]
        uwc3err = uwc3["r2se"]
        
        uwc3m = round(uwc3vals.mean(), 2)   # Mean
        
        uwc3std = round(uwc3vals.sem(), 2) * 2  # 2SE
        # UWC-1
        uwc1strings = ["uwc-1", "uwc1", "uwc_1"]
        for string in uwc1strings:        
            uwc1 = self.dataOvert.loc[self.dataOvert["sample"].str.match(string)]
            if not uwc1.empty:
                break
            
        uwc1vals = uwc1["cald11B"]
        uwc1err = uwc1["r2se"]
        uwc1m = round(uwc1vals.mean(), 2)
        uwc1std = round(uwc1vals.std() / (len(uwc1vals)**0.5), 2) * 2
        # UWC1 B/Ca = 121 mmol/mol (Foster 2013, Chem. Geo. )
        # [B] = 13 ug/g (Kasemann 2009, Chem. Geo.)
         
        # Blue calcite
        blue = self.dataOvert.loc[self.dataOvert["sample"].str.match("blue")]
        bluevals = blue["cald11B"]
        blueerr = blue["r2se"]
        bluem = round(bluevals.mean(), 2)
        bluestd = round(bluevals.std() / (len(bluevals)**0.5), 2) * 2
        
        plt.figure()
        x = range(len(uwc1vals))
        plt.scatter(x, uwc1vals, s=150, c="red", alpha=0.6, ec="k",
                    label=u"UWC-1 = {} $\pm$ {}‰".format(uwc1m, uwc1std))
        plt.errorbar(x, uwc1vals, yerr=uwc1err,
                     fmt="none", color="red", alpha=0.6)

        # Plotting uwc1 true value
        plt.hlines(7.77, -1, len(uwc1vals), colors="black")
        plt.hlines(7.77 + 0.89, -1, len(uwc1vals), colors="black", ls="--")
        plt.hlines(7.77 - 0.89, -1, len(uwc1vals), colors="black", ls="--")
        # plt.xlim(-1, len(uwcvals)+1)

        x = range(len(uwc3vals))
        plt.scatter(x, uwc3vals, s=150, c="orange", alpha=0.6, ec="k",
                    label=u"UWC-3 = {} $\pm$ {}‰".format(uwc3m, uwc3std))
        plt.errorbar(x, uwc3vals, yerr=uwc3err,
                     fmt="none", color="orange", alpha=0.6)

        plt.hlines(20.25, -1, len(uwc3vals), colors="gray")
        plt.hlines(20.25 + 0.15, -1, len(uwc3vals), colors="gray", ls="--")
        plt.hlines(20.25 - 0.15, -1, len(uwc3vals), colors="gray", ls="--")
        
        x = range(len(bluevals))
        plt.scatter(x, bluevals, s=150, c="blue", alpha=0.6, ec="k",
                    label=u"Blue = {} $\pm$ {}‰".format(bluem, bluestd))
        plt.errorbar(x, bluevals, yerr=blueerr,
                     fmt="none", color="blue", alpha=0.6)

        plt.hlines(-0.12, -1, len(bluevals), colors="gray")
        plt.hlines(-0.12 + 0.41, -1, len(bluevals), colors="gray", ls="--")
        plt.hlines(-0.12 - 0.41, -1, len(bluevals), colors="gray", ls="--")
        
        plt.ylim(-3, 25)
        plt.ylabel(r"$\delta ^{11}$B NIST951")
        plt.xlabel("Run number")
        plt.legend()
        if save:
            if name:
                plt.savefig(self.directory + "/uwc_plot_" + name + ".png")
            else:
                plt.savefig(self.directory + "/uwc_plot.png")
            
        plt.show()
        # Plot a regression between true and measured value
        plt.figure()
        plt.errorbar(-0.12, bluem, yerr=bluestd, xerr=0.41, fmt="o", mec="k")   # Blue
        plt.errorbar(7.77, uwc1m, yerr=uwc1std, xerr=0.89, fmt="o", mec="k")   # UWC1
        plt.errorbar(20.25, uwc3m, yerr=uwc3std, xerr=0.15, fmt="o", mec="k")   # UWC3
        px = np.linspace(-1, 21)
        plt.plot(px, px, "k", label="1:1 line")
        m, c = np.polyfit([-0.12, 7.77, 20.25], [bluem, uwc1m, uwc3m],1)
        x = np.linspace(-0.2, 20.5)
        plt.plot(x, x*m+c, "k:", label="measured fit")
        plt.text(1, 15, f"f(x) = {m:0.2f} * x + {c:0.2f}")
        plt.xlim(-1, 24)
        plt.ylim(-1, 24)
        plt.xlabel(r"Solution measured $\delta ^{11}$B")
        plt.ylabel(r"Laser measured $\delta ^{11}$B")
        plt.legend()
        if save:
            plt.savefig(self.directory + "/uwc_reg_plot.png")
        plt.show()
        
        plt.figure()
        plt.errorbar(blue["11/10.035"], blue["cald11B"], yerr=blue["r2se"], 
                     fmt="o", label="Blue")
        plt.errorbar(uwc1["11/10.035"], uwc1["cald11B"], yerr=uwc1["r2se"], 
                     fmt="o", label="UWC-1")
        plt.errorbar(uwc3["11/10.035"], uwc3["cald11B"], yerr=uwc3["r2se"], 
                     fmt="o", label="UWC-3") 
        plt.xlabel("11/10.035")
        plt.ylabel(r"$\delta ^{11}$B")
        plt.legend()
        plt.show()
    
    # TODO Proper [B] / B/Ca calculations 
    # TODO accuracy check of [B] / B/Ca
        
    def savedata(self, name=False):
        """Save the standards as csv."""
        self.standard_detect()
        if name:
        # Make a new dir where the rawfiles are
            ndir = self.directory + "/results_" + name + "/"
        else:
            ndir = self.directory + "/results/"
        # Try to make the directory, but skips if it already exists
        try:
            os.mkdir(ndir)  # Folder for all the reduced data
            os.mkdir(ndir + "stds")  # Folder for the individual standards
        except FileExistsError:
            pass
        # Create an array filled with "False", then add a bool array of where 
        # a certain standard is in the dataset
        stdmask = np.full(len(self.dataOvert), False)
        # Save each standards data as csv
        for st in self.stds:
            # Create a DataFrame to store and combine the data for each standard
            mask = self.dataOvert["sample"].str.match(st)
            standard_df = self.dataOvert.loc[mask]
            # Append the location mask (bool array)
            stdmask += mask
            # Save it as csv in the new directory
            standard_df.to_csv(f"{ndir}stds/{st}.csv", index=None)
        
        # Save all the standards as a single csv file
        self.dataOvert.loc[stdmask].to_csv(ndir + "standards.csv", index=None)
        # Save all the samples as a single csv file
        self.dataOvert.loc[~stdmask].to_csv(ndir + "samples.csv", index=None)
        # Save the whole sequence
        self.dataOvert.to_csv(ndir + "altogether.csv", index=None)
        # Save the data dictionary (more information)
        with open(ndir + "datadi_save.pkl", "wb") as f:
            pickle.dump(self.datadi, f)
        
        print(f"All saved in {self.directory + '/results/'} ! :) ")
        
        
        
    def save_di(self, fn):
        """Save the datadi attribute to a pickle object."""
        with open(fn, "wb") as f:
            pickle.dump(self.datadi, f)
        
        
    def doall(self):
        """Do all of the steps."""
        self.logtime()
        self.bg_sub()
        self.nistcor()
        self.ca_cor(True, True)

if __name__ == "__main__":
    # If the code is executed as a script
    from boronredc_nospot import Reduced
    
    try:  # Make sure to use QT backend for plots for manual selections
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='qt')
    except:
        pass   
    
    red = Reduced()
    red.logtime(True)
    # red.bg_sub()
    # red.nistcor()
    # red.ca_cor()
    # red.savedata()

