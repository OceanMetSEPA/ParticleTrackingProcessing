import numpy as np
import datetime
import re
import scipy.io
import warnings
from .parse_xml_row import parse_xml_row

def xml2dict(xmlFileName=None, matfileName=None, hourly_only=True):
    """ Converts an XML file into a structured dictionary format.
    
    Parameters:
        xmlFileName (str): Path to the XML file.
        matfileName (str): Optional output .mat filename.
        hourly_only (bool): If True (default), only hourly timesteps are kept.
        
    Returns:
        xmlDict (dict): Structured dictionary of extracted data.
        matfileName (str): Output .mat file name used.
    """
    
    if xmlFileName is None:
        print('Usage: xml2dict(xmlFileName)')
        return

    # Read header lines
    with open(xmlFileName, 'r') as file:
        txt = [file.readline().strip() for _ in range(100)]

    # Extract start and end times
    start_row = next((line for line in txt if 'StartTime' in line), None)
    _, start_string = parse_xml_row(start_row) if start_row else (None, None)
    end_row = next((line for line in txt if 'EndTime' in line), None)
    _, end_string = parse_xml_row(end_row) if end_row else (None, None)
    
    start_time = datetime.datetime.strptime(start_string, "%Y-%m-%d %H:%M:%S") if start_string else None
    end_time = datetime.datetime.strptime(end_string, "%Y-%m-%d %H:%M:%S") if end_string else None
    
    if start_time and end_time:
        model_duration = (end_time - start_time).total_seconds() / (24 * 3600)
    else:
        model_duration = None

    # Time step spacing from XML
    str_tag = 'TimeStep nr'
    m = [s for s in txt if str_tag in s]
    str2remove = [str_tag, '=', '<', '>', '"']
    for r in str2remove:
        m = [s.replace(r, "") for s in m]
    m = np.array([int(s) for s in m])
    timestep_offset = m.min()
    dm = np.diff(m)
    udm = np.unique(dm)
    if len(udm) != 1:
        raise ValueError('Non-equidistant timesteps!')
    else:
        dtXML = udm[0]
        print(f'Time spacing = {dtXML}, Timestep offset = {timestep_offset}')

    # Extract variable names (codes)
    codes = [re.sub(r"</?code>", "", line) for line in txt if '<code>' in line and '_' not in line]

    print(f"Start time: {start_time}")
    print(f"End time: {end_time}")
    print(f"Model duration: {model_duration} days")
    print(f"Extracted {len(codes)} codes: {codes}")

    # Always use dynamic lists now (no preallocation)
    xmlDict = {"dateTime": []}
    for code in codes:
        xmlDict[code] = []

    # Parse the full file
    with open(xmlFileName, 'r') as file:
        time_is_hourly = False
        time_index_hourly = -1
        time_index_full = -1
        particle_index = 0

        for line in file:
            var_name, val = parse_xml_row(line)

            if var_name.lower() == "timestep nr":
                raw_time_index = int(val)
                time_index_full = int((raw_time_index - timestep_offset) / dtXML)
                time_is_hourly = False

            elif var_name == "DateTime":
                dt = datetime.datetime.strptime(val, "%Y-%m-%d %H:%M:%S")
                if hourly_only:
                    if dt.minute == 0 and dt.second == 0:
                        time_is_hourly = True
                        time_index_hourly += 1

                        # ✅ Print if midnight
                        if dt.hour == 0:
                            print(f"Midnight timestep at: {dt.date()}")

                        datenum = dt.toordinal() + 366 + (dt.hour / 24)
                        xmlDict["dateTime"].append(datenum)

                        for code in codes:
                            xmlDict[code].append([])
                else:
                    datenum = (dt.toordinal() + 366 +
                               dt.hour / 24 +
                               dt.minute / 1440 +
                               dt.second / 86400)
                    # expand lists dynamically
                    while len(xmlDict["dateTime"]) <= time_index_full:
                        xmlDict["dateTime"].append(np.nan)
                        for code in codes:
                            xmlDict[code].append([])
                    xmlDict["dateTime"][time_index_full] = datenum

            elif var_name == "Particle Nr":
                particle_index = int(val) - 1

            elif var_name in codes:
                if hourly_only and time_is_hourly:
                    while len(xmlDict[var_name][time_index_hourly]) <= particle_index:
                        xmlDict[var_name][time_index_hourly].append(np.nan)
                    xmlDict[var_name][time_index_hourly][particle_index] = float(val)
                elif not hourly_only:
                    while len(xmlDict[var_name][time_index_full]) <= particle_index:
                        xmlDict[var_name][time_index_full].append(np.nan)
                    xmlDict[var_name][time_index_full][particle_index] = float(val)

    # Finalize: convert lists → arrays
    xmlDict["dateTime"] = np.array(xmlDict["dateTime"]).reshape(-1, 1)
    for code in codes:
        # ragged to padded
        maxNp = max(len(p) for p in xmlDict[code])
        Nt = len(xmlDict[code])
        arr = np.full((maxNp, Nt), np.nan)
        for t, plist in enumerate(xmlDict[code]):
            for p, val in enumerate(plist):
                arr[p, t] = val
        xmlDict[code] = arr

    print("Saving matfile...")
    try:
        if matfileName is None:
            matfileName = xmlFileName.replace('.xml', '.mat')
        scipy.io.savemat(matfileName, xmlDict)
    except Exception as e:
        warnings.warn(f"Something went wrong creating .mat: {e}")

    print("Done.")
    return xmlDict, matfileName
