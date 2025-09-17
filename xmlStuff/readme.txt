Python code for converting .xml file to .matfile

NB, code developed to convert xml file with single source/class. Will need to be modified if multiple sources included in file.

To use:

1. Open Spyder.
2. Make sure the xmltools folder is somewhere in Python’s path (e.g., C:\Python\).
3. Open run_xml_processing.py in Spyder.
4. Change model_path if needed.
5. Press Run File (F5).

This will process all subfolders containing an xml file (ignoring those that already have an associated matfile).

Other functions:
process_xml_folder - process a single folder containing xml
xml2dict - actually does the conversion from .xml to python dict
parse_xml_row - called by function above; extracts values from formatted row in xml file. Shouldn't need to be called directly. 

