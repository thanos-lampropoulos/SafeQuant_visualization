import streamlit as st
import pandas as pd
import plotly.express as px
import re
import math

###############
# defining a function

def sq_processing(sq_df):

    """
    This function processes the PROTEIN.tsv file that SafeQuant serves as output and returns a dictionary whose values are Pandas dataframes.

    Parameters
    ----------
    Safequant output (PROTEIN.tsv)
        The tsv file that SafeQuant returns with the protein data.
    
    Returns
    -------
    dictionary : dict
        The collection of Pandas dataframes.
    """

    df = sq_df.copy()
    
    # creating the list of column names
    # creating the list of column names to keep
    # creating the list of column names to drop
    
    # the * operator unpacks the .columns set, the list() method would also work 
    column_names = [*df.columns]
    
    # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
    columns_to_keep = [name for name in column_names if re.search( r"(proteinName)|(^ac)|(geneName)|(proteinDescription)|(nbPeptides)|(^pValue)|(^qValue)|(^log2ratio)", name)]
    
    # the columns_to_drop list gets processed and complete with the for loop below
    columns_to_drop = column_names.copy()
    
    for column in columns_to_keep:
    
        columns_to_drop.remove(column)

    # dropping the unnecessary columns

    df.drop(labels = columns_to_drop, axis = 1, inplace = True)

    # renaming columns

    df.rename(columns = {"proteinName" : "Protein Name",
                              "ac" : "Accession",
                              "geneName" : "Gene Name",
                              "proteinDescription" : "Protein Description"
                             }, inplace = True)

    # shortening the info in the "Protein Description" column

    df["Protein Description"] = df["Protein Description"].apply(lambda x: re.sub(r"\sOS=.+$", "", x))
    
    # removing additional accessions from the "Protein Name" column
    
    df["Protein Name"] = df["Protein Name"].apply(lambda x: re.sub(r";.+$", "", x))
    
    # creating a new column with a shortened "Protein Name" and making it the second column of the data frame
    # ATTENTION: this will disguise ligands coming from another species!!!
    
    df["Protein Name (short)"] = df["Protein Name"].apply(lambda x: re.sub(r"(^sp\|.+\|)|(_.+$)", "", x))
    
    column_to_move = df.pop("Protein Name (short)")
    
    df.insert(1, "Protein Name (short)", column_to_move)

    # selecting the columns that will be logarithmized

    # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
    columns_for_log = [name for name in columns_to_keep if re.search( r"(^qValue)", name)]
    
    # logarithmizing the qValues
    for name in columns_for_log:
        
        df[f"-log10({name})"] = df[f"{name}"].apply(lambda x: abs(math.log10(x)))

    # selecting columns whose names end with all possible ligands, e.g. the log2 columns
    # the columns_for_log list from above can also be used instead

    columns_ligand = [name for name in df.columns if re.search(r"^log2", name)]

    # finding out what are the names of the treatment arms
    # IMPORTANT: the name of the treatment arm must not contain "_" and always be in the format e.g. log2ratio_NAME
    
    treatment_arms =[]
    
    for name in columns_ligand:
    
        # this returns the match object for each treatment arm
        match = re.search(r"_.+$", name)
    
        # this returns the matching string for each treatment arm after stripping the "_"
        treatment_arms.append(match.group().lstrip("_"))

    # these columns need always be present
    columns_obligatory = [*df.iloc[:, 0:6].columns]
    
    # this placeholder dictionary will store each ligand's dataframe
    df_collection = {}
    
    for arm in treatment_arms:
    
        # create a dataframe for each treatment arm with the columns to be kept in each iteration, i.e. obligatory columns, plus ligand-specific columns
        # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
        columns_to_include = [name for name in df.columns if name in columns_obligatory or re.search(fr"_{arm}\)?$", name)]
        
        df_interim = df.loc[:, columns_to_include].copy()
        
        df_collection[f"{arm}"] = df_interim.copy()
    
        ################
        # export the dataframes as tsv files
        # renaming the -log10(qValue) column before exporting for compatibility with Excel
    
        #df_interim.rename(columns = {f"-log10(qValue_{arm})" : f"'-log10(qValue_{arm})"}, inplace = True)
    
        #df_interim.to_csv(f'{ligand_on_the_left}_vs_{arm}_{number_of_peptides}.tsv', sep = '\t', index=False)
        
        ################
        
    # this dictionary contains the final tables for each ligand/treatment arm 
    return df_collection  
    

###############
# defining a function

def sq_processing_manual(sq_df):

    """
    This function processes the PROTEIN.tsv file that SafeQuant serves as output and returns a dictionary whose values are Pandas dataframes.
    This is a copy of function sq_processing(), but also returns the Pandas dataframes as tsv files.
    The tsv files are meant to be used with Excel for manual visualization.
    
    Parameters
    ----------
    Safequant output (PROTEIN.tsv)
        The tsv file that SafeQuant returns with the protein data.
    
    Returns
    -------
    dictionary : dict
        The collection of Pandas dataframes.
    tsv files: tsv
        Each Pandas dataframe is exported as a tsv file.
    """

    df = sq_df.copy()
    
    # creating the list of column names
    # creating the list of column names to keep
    # creating the list of column names to drop
    
    # the * operator unpacks the .columns set, the list() method would also work 
    column_names = [*df.columns]
    
    # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
    columns_to_keep = [name for name in column_names if re.search( r"(proteinName)|(^ac)|(geneName)|(proteinDescription)|(nbPeptides)|(^pValue)|(^qValue)|(^log2ratio)", name)]
    
    # the columns_to_drop list gets processed and complete with the for loop below
    columns_to_drop = column_names.copy()
    
    for column in columns_to_keep:
    
        columns_to_drop.remove(column)

    # dropping the unnecessary columns

    df.drop(labels = columns_to_drop, axis = 1, inplace = True)

    # renaming columns

    df.rename(columns = {"proteinName" : "Protein Name",
                              "ac" : "Accession",
                              "geneName" : "Gene Name",
                              "proteinDescription" : "Protein Description"
                             }, inplace = True)

    # shortening the info in the "Protein Description" column

    df["Protein Description"] = df["Protein Description"].apply(lambda x: re.sub(r"\sOS=.+$", "", x))
    
    # removing additional accessions from the "Protein Name" column
    
    df["Protein Name"] = df["Protein Name"].apply(lambda x: re.sub(r";.+$", "", x))
    
    # creating a new column with a shortened "Protein Name" and making it the second column of the data frame
    # ATTENTION: this will disguise ligands coming from another species!!!
    
    df["Protein Name (short)"] = df["Protein Name"].apply(lambda x: re.sub(r"(^sp\|.+\|)|(_.+$)", "", x))
    
    column_to_move = df.pop("Protein Name (short)")
    
    df.insert(1, "Protein Name (short)", column_to_move)

    # selecting the columns that will be logarithmized

    # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
    columns_for_log = [name for name in columns_to_keep if re.search( r"(^qValue)", name)]
    
    # logarithmizing the qValues
    for name in columns_for_log:
        
        df[f"-log10({name})"] = df[f"{name}"].apply(lambda x: abs(math.log10(x)))

    # selecting columns whose names end with all possible ligands, e.g. the log2 columns
    # the columns_for_log list from above can also be used instead

    columns_ligand = [name for name in df.columns if re.search(r"^log2", name)]

    # finding out what are the names of the treatment arms
    # IMPORTANT: the name of the treatment arm must not contain "_" and always be in the format e.g. log2ratio_NAME
    
    treatment_arms =[]
    
    for name in columns_ligand:
    
        # this returns the match object for each treatment arm
        match = re.search(r"_.+$", name)
    
        # this returns the matching string for each treatment arm after stripping the "_"
        treatment_arms.append(match.group().lstrip("_"))

    # these columns need always be present
    columns_obligatory = [*df.iloc[:, 0:6].columns]
    
    # this placeholder dictionary will store each ligand's dataframe
    df_collection = {}
    
    for arm in treatment_arms:
    
        # create a dataframe for each treatment arm with the columns to be kept in each iteration, i.e. obligatory columns, plus ligand-specific columns
        # Match objects always have a boolean value of True. re.search() returns Match objects if it finds a match
        columns_to_include = [name for name in df.columns if name in columns_obligatory or re.search(fr"_{arm}\)?$", name)]
        
        df_interim = df.loc[:, columns_to_include].copy()
        
        df_collection[f"{arm}"] = df_interim.copy()
    
        ################
        # export the dataframes as tsv files
        # renaming the -log10(qValue) column before exporting for compatibility with Excel
    
        df_interim.rename(columns = {f"-log10(qValue_{arm})" : f"'-log10(qValue_{arm})"}, inplace = True)
    
        df_interim.to_csv(f'{ligand_on_the_left}_vs_{arm}_{number_of_peptides}.tsv', sep = '\t', index=False)
           
        with open(f'{ligand_on_the_left}_vs_{arm}_{number_of_peptides}.tsv', mode = 'rb') as f:
            st.download_button(label=f'Download {ligand_on_the_left}_vs_{arm}_{number_of_peptides}.tsv', data = f, file_name = f'{ligand_on_the_left}_vs_{arm}_{number_of_peptides}.tsv', mime= 'application/octet-stream')
        
        
        ################
        
    # this dictionary contains the final tables for each ligand/treatment arm 
    return df_collection    
    
    
###############
# defining a function

def sq_plot(dictionary, enrichment_threshold, statistical_threshold):
    
    """
    This function visualizes the elements of a dictionary whose values are Pandas dataframes containing data from SafeQuant.

    Parameters
    ----------
    dictionary : dict
        The collection of Pandas dataframes.
    
    Returns
    -------
    Plotly Plots
        The plots created by Plotly.
    """
   
    # reminder: for loops with dictionaries in python loop through the keys
    # the key needs to be used as dictionary[key] in the for loop in order to get the value
    
    for key in dictionary:
        
        ################################################## 
        # determine the range of the x axis based on log2 column
        
        log2max = dictionary[key].iloc[:, 6].max()

        log2min = abs(dictionary[key].iloc[:, 6].min())
        

        if log2max > log2min:

            log2range = math.ceil(log2max + 1)

            if log2range % 2 == 0:
                
                pass
          
            else:

                log2range = log2range + 1

        elif log2max < log2min:

            log2range = math.ceil(log2min + 1)

            if log2range % 2 == 0:

                pass
          
            else:

                log2range = log2range + 1
                
        # this is for the rare case were log2max and log2min might be equal
        
        else:

            log2range = math.ceil(log2min + 1)

            if log2range % 2 == 0:

                pass
          
            else:

                log2range = log2range + 1

        ##################################################
        # determine the range of the y axis based on the -log10 column
        
        log10max = dictionary[key].iloc[:, 9].max()

        log10range = math.ceil(log10max + 2)

        ##################################################
        # using plotly to draw the volcano plot for each pairwise comparison
        
        fig = px.scatter(dictionary[key],
                         x = dictionary[key].iloc[: , 6],
                         y = dictionary[key].iloc[: , 9], 
                         hover_name = dictionary[key].iloc[: , 0],
                         hover_data = [dictionary[key].iloc[: , 5], dictionary[key].iloc[: , 6], dictionary[key].iloc[: , 9]],
                         labels = {dictionary[key].iloc[: , 6].name : "log\u2082(fold change)", dictionary[key].iloc[: , 9].name : "-log\u2081\u2080(adjusted p-value)"}
                         #text = dictionary[key].iloc[: , 1]
                        )
                 
        fig.update_traces(textposition='top center')

        fig.update_traces(marker={"size" : 6,
                                  "line": {"width" : 1, "color" : "black"},
                                  "color" : "teal"
                                 }
                         )

        # setting background properties
        fig.update_layout(plot_bgcolor = "white")

        # setting title properties
        fig.update_layout(title_text = f"{ligand_on_the_left} vs {key} ({project_info}, {number_of_peptides})")
        fig.update_layout(title = {'x' : 0.5, 'y' : 0.96,'xanchor' : 'center', 'yanchor' : 'top'})        
        fig.update_xaxes(title_font = {"size": 16},  title_standoff = 10)
        fig.update_yaxes(title_font = {"size": 16},  title_standoff = 10)
        
        # setting axes range and tick properties
        fig.update_xaxes(range = [-log2range, log2range], fixedrange = False, dtick = 2, ticklabelstandoff = 7)
        fig.update_yaxes(range = [0, log10range], fixedrange = False, dtick = 1, ticklabelstandoff = 7)
        fig.update_yaxes(ticklabelstep=1)

        # setting axes line properties
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror = True)
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror = True)
        
        # setting zero line properties
        fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='black')
        fig.update_yaxes(zeroline=False, zerolinewidth=3, zerolinecolor='black')

        # setting grid properties (Plotly accepts hex colours as strings)
        fig.update_yaxes(showgrid = True, gridcolor = '#bbbbbf', gridwidth = 1)
        fig.update_xaxes(showgrid = True, gridcolor = '#bbbbbf', gridwidth = 1)

        # adding background shapes
        fig.add_shape(type = "rect",
                      x0 = enrichment_threshold, 
                      y0 = 0, 
                      x1 = log2range, 
                      y1 = statistical_threshold,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor = "blue",
                      opacity = 0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = 0, 
                      y0 = 0, 
                      x1 = enrichment_threshold, 
                      y1 = log10range,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor = "blue",
                      opacity=0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = -enrichment_threshold, 
                      y0 = 0, 
                      x1 = -log2range, 
                      y1 = statistical_threshold,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor="blue",
                      opacity=0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = 0, 
                      y0 = 0, 
                      x1 = -enrichment_threshold, 
                      y1 = log10range,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor="blue",
                      opacity=0.1
                     )
               
        #fig.show()
        fig.write_html(f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}.html", auto_open=False)
        st.plotly_chart(fig, theme = None)
        
        with open(f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}.html", mode = 'rb') as f:
            st.download_button(label=f"Download {project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}.html", data = f, file_name = f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}.html", mime= 'application/octet-stream')

###############
# defining a function

def sq_plot_text(dictionary, enrichment_threshold, statistical_threshold):
    
    """
    This function visualizes the elements of a dictionary whose values are Pandas dataframes containing data from SafeQuant.
    This is a copy of function sq_plot(), but shows by default text annotations in the volcano plots.

    Parameters
    ----------
    dictionary : dict
        The collection of Pandas dataframes.
    
    Returns
    -------
    Plotly Plots
        The plots created by Plotly.
    """
   
    # reminder: for loops with dictionaries in python loop through the keys
    # the key needs to be used as dictionary[key] in the for loop in order to get the value
    
    for key in dictionary:
        
        ################################################## 
        # determine the range of the x axis based on log2 column
        
        log2max = dictionary[key].iloc[:, 6].max()

        log2min = abs(dictionary[key].iloc[:, 6].min())
        

        if log2max > log2min:

            log2range = math.ceil(log2max + 1)

            if log2range % 2 == 0:
                
                pass
          
            else:

                log2range = log2range + 1

        elif log2max < log2min:

            log2range = math.ceil(log2min + 1)

            if log2range % 2 == 0:

                pass
          
            else:

                log2range = log2range + 1
                
        # this is for the rare case were log2max and log2min might be equal
        
        else:

            log2range = math.ceil(log2min + 1)

            if log2range % 2 == 0:

                pass
          
            else:

                log2range = log2range + 1

        ##################################################
        # determine the range of the y axis based on the -log10 column
        
        log10max = dictionary[key].iloc[:, 9].max()

        log10range = math.ceil(log10max + 2)

        ##################################################
        # using plotly to draw the volcano plot for each pairwise comparison
        
        fig = px.scatter(dictionary[key],
                         x = dictionary[key].iloc[: , 6],
                         y = dictionary[key].iloc[: , 9], 
                         hover_name = dictionary[key].iloc[: , 0],
                         hover_data = [dictionary[key].iloc[: , 5], dictionary[key].iloc[: , 6], dictionary[key].iloc[: , 9]],
                         labels = {dictionary[key].iloc[: , 6].name : "log\u2082(fold change)", dictionary[key].iloc[: , 9].name : "-log\u2081\u2080(adjusted p-value)"},
                         text = dictionary[key].iloc[: , 1]
                        )
                 
        fig.update_traces(textposition='top center')

        fig.update_traces(marker={"size" : 6,
                                  "line": {"width" : 1, "color" : "black"},
                                  "color" : "teal"
                                 }
                         )

        # setting background properties
        fig.update_layout(plot_bgcolor = "white")

        # setting title properties
        fig.update_layout(title_text = f"{ligand_on_the_left} vs {key} ({project_info}, {number_of_peptides})")
        fig.update_layout(title = {'x' : 0.5, 'y' : 0.96,'xanchor' : 'center', 'yanchor' : 'top'})        
        fig.update_xaxes(title_font = {"size": 16},  title_standoff = 10)
        fig.update_yaxes(title_font = {"size": 16},  title_standoff = 10)
        
        # setting axes range and tick properties
        fig.update_xaxes(range = [-log2range, log2range], fixedrange = False, dtick = 2, ticklabelstandoff = 7)
        fig.update_yaxes(range = [0, log10range], fixedrange = False, dtick = 1, ticklabelstandoff = 7)
        fig.update_yaxes(ticklabelstep=1)

        # setting axes line properties
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror = True)
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror = True)
        
        # setting zero line properties
        fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='black')
        fig.update_yaxes(zeroline=False, zerolinewidth=3, zerolinecolor='black')

        # setting grid properties (Plotly accepts hex colours as strings)
        fig.update_yaxes(showgrid = True, gridcolor = '#bbbbbf', gridwidth = 1)
        fig.update_xaxes(showgrid = True, gridcolor = '#bbbbbf', gridwidth = 1)

        # adding background shapes
        fig.add_shape(type = "rect",
                      x0 = enrichment_threshold, 
                      y0 = 0, 
                      x1 = log2range, 
                      y1 = statistical_threshold,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor = "blue",
                      opacity = 0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = 0, 
                      y0 = 0, 
                      x1 = enrichment_threshold, 
                      y1 = log10range,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor = "blue",
                      opacity=0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = -enrichment_threshold, 
                      y0 = 0, 
                      x1 = -log2range, 
                      y1 = statistical_threshold,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor="blue",
                      opacity=0.1
                     )
        
        fig.add_shape(type="rect",
                      x0 = 0, 
                      y0 = 0, 
                      x1 = -enrichment_threshold, 
                      y1 = log10range,
                      line = {"color" : "blue", "width" : 0},
                      fillcolor="blue",
                      opacity=0.1
                     )
               
        #fig.show()
        fig.write_html(f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}_withText.html", auto_open=False)
        st.plotly_chart(fig, theme = None)        
        
        with open(f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}_withText.html", mode = 'rb') as f:
            st.download_button(label=f"Download {project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}_withText.html", data = f, file_name = f"{project_info}_{number_of_peptides}_{ligand_on_the_left}_vs_{key}_withText.html", mime= 'application/octet-stream')


########################################
# Introduction to the app.
st.title("Processing and visualizing SafeQuant results.")

st.write("--------------------------------------------------")

st.write("""Enter the details of your project and then upload the PROTEIN.tsv file from Safequant.

After uploading a PROTEIN.tsv file:

a. the results will be processed and you can download a tsv file for each pairwise comparison.

b. an interactive volcano plot (without text annotations) will be created, which can be downloaded.

c. an interactive volcano plot (with text annotations) will be created, which can be downloaded.""")

st.write("--------------------------------------------------")

########################################
# Defining project-wide variables (just strings).
# These variables need always be defined manually.

project_info = st.text_input(
    label = "Enter here the name of the project (e.g. R399).",
    value = "",
    key = "project_info_key") #e.g. "R424"

if project_info:
    st.write(f"The name of the project is: **{st.session_state["project_info_key"]}**.")

ligand_on_the_left = st.text_input(
    label = "Enter here the name of the ligand that will be depicted on the left of the volcano plot (e.g. TRFE).",
    value = "",
    key = "ligand_on_the_left_key") # e.g."IL38"

if ligand_on_the_left:
    st.write(f"The name of the ligand on the left is: **{st.session_state["ligand_on_the_left_key"]}**.")

number_of_peptides = st.text_input(
    label = "The number of peptides (input as 1pep or 2pep).",
    value = "",
    key = "number_of_peptides_key") #e.g. "2pep"

if number_of_peptides:
    st.write(f"The number of peptides is: **{st.session_state["number_of_peptides_key"]}**.")

st.write("--------------------------------------------------")

########################################
# Uploading the PROTEIN.tsv file from SafeQuant.

# this session_state key serves as a placeholder for "file" below
# this key will be changed every time the "reset" button is clicked (see below),
# thereby ensuring that the "file" widget gets replaced by a new empty instance
if "file_uploader_key" not in st.session_state:
    st.session_state["file_uploader_key"] = 0

file = st.file_uploader(
    label = "Select the PROTEIN.tsv file.", 
    type = "tsv",
    key = st.session_state["file_uploader_key"])

if file is not None:

    dual = pd.read_csv(file, sep = '\t')

    dual_df = dual.copy()

    df_for_review = dual_df.head(n = 5).copy()
    
    st.write("This is how the data you uploaded looks like:")
        
    st.dataframe(df_for_review)
        
    st.write("--------------------------------------------------")

    ########################################

    # processing alternative 1: processing the SafeQuant tsv file

    #dict_for_viz = sq_processing(dual_df)


    # alternative 2: processing the SafeQuant tsv file and saving the tsv files
    
    
    st.write("#### You can download the processed results as tsv files.")
    
    
    dict_for_viz = sq_processing_manual(dual_df)
    
    
    st.write("--------------------------------------------------")
    
    # setting the enrichment and statistical significance thresholds
    
    st.write("#### You can set the enrichment and statistical significance thresholds.")
    
    enrichment_thr = st.slider(
        label = "Set the enrichment threshold (log2 space):",
        min_value = 0.0,
        max_value = 10.0,
        value = 2.0,
        step = 0.5)
    
    statistical_thr = st.slider(
        label = "Set the statistical threshold (-log10 space):",
        min_value = 0.0,
        max_value = 3.0,
        value = 2.0,
        step = 0.1)

    
    # visualization alternative 1: plotly plots without text annotations
    
    
    st.write("#### You can download the volcano plots without annotations.")
    
    
    sq_plot(dict_for_viz, enrichment_thr, statistical_thr)
    
    
    st.write("--------------------------------------------------")

    
    # visualization alternative 2: plotly plots with text annotations
    
    
    st.write("#### You can download the volcano plots with annotations.")
    
    
    sq_plot_text(dict_for_viz, enrichment_thr, statistical_thr)
    
    
    st.write("--------------------------------------------------")
    st.write("--------------------------------------------------")
    
else:
    st.write("### Please upload a file to use the app.")


# creating a button to reset the app

def reset():
    
    st.session_state["number_of_peptides_key"] = ""
     
    st.session_state["ligand_on_the_left_key"] = ""
    
    st.session_state["project_info_key"] = ""
    
    st.session_state["file_uploader_key"] += 1

st.button(
    "Click here to reset the app (or reload the webpage instead).",
    on_click = reset)

#st.session_state

