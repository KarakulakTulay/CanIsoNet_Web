from flask import Flask, render_template, request
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import json
from flaskext.mysql import MySQL
import statistics
import yaml
from flask_paginate import Pagination, get_page_parameter
#from flask_paginate import get_page_parameter

app = Flask(__name__)

db = yaml.load(open('/home/abxka/CanIsoNet_Web/db.yaml'), Loader=yaml.FullLoader)
app.config['MYSQL_DATABASE_HOST'] = db['mysql_host']
app.config['MYSQL_DATABASE_USER'] = db['mysql_user']
app.config['MYSQL_DATABASE_PASSWORD'] = db['mysql_password']
app.config['MYSQL_DATABASE_DB'] = db['mysql_db']

mysql = MySQL()
mysql.init_app(app)

@app.route("/")
def home():

    # get cancer types
    cur = mysql.get_db().cursor()
    cur.execute('''SELECT * FROM cancertypesinproject ORDER BY `Disease Name`;''')
    cancer_types = cur.fetchall()
    cancertypes = pd.DataFrame(cancer_types,columns=['PCAWG-Code', 'Cancer Type Name', 'PCAWG_GTEx', 'Dataset'])
    cancertypes_dict = cancertypes.to_dict(orient='records')

    # get available human transcripts in the dataset
    cur2 = mysql.get_db().cursor()
    cur2.execute(''' SELECT DomCancerTrans, GeneName1_x, ENSG, Transcript_Name FROM ENST_Genename_ENSG_TranscriptName ''')
    genename_cmdt = cur2.fetchall()
    genenamecmdt = pd.DataFrame(genename_cmdt, columns=['DomCancerTrans','GeneName1_x', 'ENSG', 'Transcript_Name'])
    genenamecmdt_gene_list = list(genenamecmdt.iloc[:,1].unique())
    genenamecmdt_gene_cmdt = list(genenamecmdt.iloc[:,0].unique())
    genenamecmdt_gene_tn = list(genenamecmdt.iloc[:,3].unique())
    temp_dict_human = genenamecmdt.to_dict(orient='records')

    # get available mouse transcripts in the dataset
    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT * FROM ENST_Genename_ENSG_TranscriptName_OIH_only_dMDTs ''')
    genename_cmdtmouse = cur3.fetchall()
    genename_cmdt_mouse = pd.DataFrame(genename_cmdtmouse, columns=['DomCancerTrans','GeneName1_x', 'ENSG', 'Transcript_Name'])
    genename_cmdt_mouse_gene_list = list(genename_cmdt_mouse.iloc[:,1].unique())
    genename_cmdt_mouse_gene_cmdt = list(genename_cmdt_mouse.iloc[:,0].unique())
    genename_cmdt_mouse_gene_tn = list(genename_cmdt_mouse.iloc[:,3].unique())
    temp_dict_mouse = genename_cmdt_mouse.to_dict(orient='records')

    isoform_list =  genenamecmdt_gene_cmdt + genename_cmdt_mouse_gene_cmdt + genenamecmdt_gene_tn + genename_cmdt_mouse_gene_tn
    gene_list = genenamecmdt_gene_list + genename_cmdt_mouse_gene_list
    transcript_list = genenamecmdt_gene_tn + genename_cmdt_mouse_gene_tn

    return render_template('home.html', data=cancertypes_dict,  data2=gene_list, data3=isoform_list, temp_dict_human=temp_dict_human, temp_dict_mouse=temp_dict_mouse, genenamecmdt_gene_tn=transcript_list)


@app.route("/specificdMDTs", methods=['GET', 'POST'])
def specificdMDTs():

    cur2 = mysql.get_db().cursor()
    page = request.args.get(get_page_parameter(), type=int, default=1)
    row_to_display = 15
    offset = (page - 1) * row_to_display



    cur2.execute(''' SELECT DISTINCT Tissue, ENSG, GeneName1, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1,
                     DomCancerTrans FROM interactionDisruptionInDominantTranscripts_int_anno_Human UNION SELECT DISTINCT Tissue, ENSG, GeneName1, TotalNumberOfStringInt, cMDTnumberOfUniqMissedInt, Pfam1, Domain1,
                     cMDT  FROM interactionDisruptionInDominantTranscripts_Human_Disease ORDER BY Tissue ASC LIMIT %s, %s''' , (offset, row_to_display))

    total=2034*15

    isonet_tuple = cur2.fetchall()

    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1',
                                            'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans',
                                            'Pfam1', 'Domain1', 'DomCancerTrans'])


    df = df.rename(columns={'ENSG': 'ENSGid', 'NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions','TotalNumberOfStringInt': 'NumberOfStringInt'})
    df = df.drop_duplicates()
    df['Domain1'].replace({"-":"None"},inplace=True)
    df = df.rename(columns={'Tissue':'CancerType2'})
    df['Pfam1'].replace({"-":"None"},inplace=True)
    df['MissedInteractions'].replace({-1:0},inplace=True)
    df['NumberOfStringInt'].replace({-1:0},inplace=True)
    #df['NumberOfStringInt'].replace({0:1},inplace=True)

    df['Percentage'] = df['MissedInteractions'].divide(df['NumberOfStringInt'])*100
    df['Percentage'].replace(np.nan, 0, inplace=True)
    df['Percentage'].replace(-0, 0, inplace=True)
    df['Percentage'].replace(-np.nan, 0, inplace=True)
    df.replace([np.inf, -np.inf] , 0, inplace=True)
    df['Percentage'] = df['Percentage'].round(2)

    df = df.drop_duplicates()

    pagination = Pagination(page=page, per_page=row_to_display, total=total, record_name='df', bs_version=4)
    data_dict = df.to_dict(orient='records')

    return render_template('specificdMDTs.html', temp_dict2=data_dict, pagination=pagination)



@app.route("/dtMDTs")
def dtMDTs():

    # get cancer specific MDTs and their frequencies of occurance
    cur = mysql.get_db().cursor()
    cur.execute(''' SELECT CancerType, cMDT, GeneName, Count, Total, Frequency, mart_export.`Associated Transcript Name` FROM dtMDT_Frequency LEFT JOIN mart_export ON dtMDT_Frequency.cMDT = mart_export.`Ensembl Transcript ID` ORDER BY Frequency DESC''')
    TableS9_tuple = cur.fetchall()
    result = pd.DataFrame(TableS9_tuple, columns=['CancerType', 'cMDT', 'GeneName','Count',  'Total', 'Frequency', 'Transcript_Name'])
    result[['Splitted','CancerType2']] = result.CancerType.str.split('.', expand=True)
    result2 = result[['CancerType2', 'cMDT', 'GeneName','Count',  'Total', 'Frequency', 'Transcript_Name']]
    result2 = result2.astype({'Count': int, 'Total': int})

    temp_dict = result2.to_dict(orient='record')

    return render_template('dtMDTs.html', temp_dict=temp_dict)

@app.route("/download", methods=['GET', 'POST'])
def download():

    return render_template("download.html")

@app.route("/contribute", methods=['GET', 'POST'])
def contribute():

    return render_template("contribute.html")

@app.route("/help", methods=['GET', 'POST'])
def help():

    cur = mysql.get_db().cursor()
    cur.execute(''' SELECT CancerType, Total FROM sample_sizes ''')
    TableS4_tuple = cur.fetchall()
    TableS4 = pd.DataFrame(TableS4_tuple, columns=['CancerType','Total'])
    TableS4 =  TableS4.rename(columns={'CancerType': 'CancerType','Total':'Total'})
    TableS4_uniq = TableS4.drop_duplicates()

    graphJSON =  sample_size(TableS4_uniq)

    return render_template("help.html", TableS4=TableS4, graphJSON=graphJSON)

def sample_size(TableS4):

    # plot the sample size for each disease
    x = TableS4.sort_values(by=['Total']).iloc[:,0].str.split(".", n=1, expand=True).iloc[:,1]
    y = TableS4.sort_values(by=['Total']).iloc[:,1]
    data = [go.Bar(x=x, y=y, marker_color = 'indianred', name='Number of Disease Samples')]
    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route("/Cancer", methods=['GET', 'POST'])
def Cancer():

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, GeneName1, CancerSampleId, DomCancerTrans, RelMDIexpressionInCancer, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_int_anno_Human
    LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_int_anno_Human.DomCancerTrans = mart_export.`Ensembl Transcript ID`
    WHERE interactionDisruptionInDominantTranscripts_int_anno_Human.Tissue = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'RelMDIexpressionInCancer', 'Associated Transcript Name'])
    df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
    df_iso = df[['CancerType2', 'Tissue', 'CancerSampleId', 'GeneName1', 'DomCancerTrans', 'Associated Transcript Name', 'RelMDIexpressionInCancer']]
    df_iso2 = df_iso.drop_duplicates()
    result =  df_iso2.rename(columns={'Associated Transcript Name':'Transcript_Name', 'RelMDIexpressionInCancer': 'dMDTenrichment'})
    temp_dict = result.to_dict(orient='records')

    ## second table in the page representing second graph - number of MDTs in each sample
    cur2 = mysql.get_db().cursor()
    sql2 = '''SELECT CancerType, SampleID, NumberOfMDTs FROM mdts_and_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    dataset_muts = cur2.fetchall()
    muts_df = pd.DataFrame(dataset_muts, columns= ['CancerType', 'SampleID', 'NumberOfMDTs' ])
    muts_df[['Splitted','CancerType2']] = muts_df.CancerType.str.split('.', expand=True)
    dataset = muts_df.to_dict(orient='records')

    return render_template("Cancer_Based.html", SampleCancerType = SampleCancerType, data = temp_dict, data2 = dataset)

def CancerSpecific(SampleCancerType):
    # take the cancer type from html
    SampleCancerType = request.args.get('cancer')

    # plot the first 10 mostly occuring dMDTs in the cancer type
    cur = mysql.get_db().cursor()
    sql = ''' SELECT cMDT, Frequency, CancerType, mart_export.`Associated Transcript Name` FROM tables4 LEFT JOIN mart_export ON tables4.cMDT = mart_export.`Ensembl Transcript ID`
            WHERE tables4.CancerType = %s ORDER BY Frequency DESC LIMIT 10 '''
    adr = (SampleCancerType, )
    cur.execute(sql,adr)
    TableS4_tuple = cur.fetchall()
    result = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','CancerType', 'Transcript_Name'])

    cur2 = mysql.get_db().cursor()
    sql2 = ''' SELECT CancerType, SampleID, NumberOfMDTs  FROM mdts_and_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    cMDT_mdts_muts = cur2.fetchall()
    cMDT_dist = pd.DataFrame(cMDT_mdts_muts, columns=['CancerType', 'SampleID','NumberOfMDTs'])



    # extract dMDT TPMs values from pcawg data - Top 10 dMDT
    enst_list = result['cMDT'].values.tolist()
    tissue = request.args.get('cancer')
    pcawg = '_pcawg'

    tissuetype = tissue.split('.')[1].replace('-','_')

    enst_tpms = mysql.get_db().cursor()
    my_tissue2 = tissuetype + pcawg
    result_enst_tpms = pd.DataFrame()
    for eachdMDT in enst_list:
        enst_tpms_sql = ''' SELECT * FROM ''' + my_tissue2 + ''' WHERE Feature REGEXP %s '''
        enst_tpms_adr = (eachdMDT,)
        enst_tpms.execute(enst_tpms_sql, enst_tpms_adr)
        enst_tpms_tuple = enst_tpms.fetchall()
        result_enst_tpms_df = pd.DataFrame(enst_tpms_tuple)
        result_enst_tpms = result_enst_tpms.append(result_enst_tpms_df, ignore_index = True)

    medians_each_row_list = []
    result_enst_tpms_TPMs = result_enst_tpms.iloc[:,1:]

    medians_each_row = result_enst_tpms_TPMs.median(axis=1)

    medians_each_row_list = medians_each_row.values.tolist()
    medians_rounds = [round(number, 2) for number in medians_each_row_list]
    medians_rounds_text = ["Median TPM: " + str(item) for item in medians_rounds]

    plot = make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3], subplot_titles=("Frequency of top 10 dMDT with median TPM values", "Number of dMDT"))


    trace1 = go.Bar(
                #name = '% of ENSTs across ' + SampleCancerType,
                marker_color = 'indianred',
                x=result.iloc[:,3],
                y=result.iloc[:,1]*100,
                text = medians_rounds,
                hovertext = medians_rounds_text,
                textposition='inside'
            )

    trace2 = go.Box(y=cMDT_dist.NumberOfMDTs,boxpoints='outliers',jitter=0.3,
                    marker_color = '#00CC96', name=' ')

    plot.append_trace(trace1, row=1, col=1)
    plot.update_yaxes(title_text="Frequency of Transcripts in Samples (%)", row=1, col=1)
    plot.append_trace(trace2, row=1, col=2)
    plot.update_yaxes(title_text="Distribiton of dMDT Across Samples", row=1, col=2)
    plot.update_xaxes(title_text=' ', row=1, col=2)
    plot.update_layout(showlegend=False)

    graphJSON = json.dumps(plot, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/CancerSpecific', methods=['GET', 'POST'])
def change_features5():

    SampleCancerType = request.args.get('cancer')
    graphJSON = CancerSpecific(SampleCancerType)

    return graphJSON

@app.route("/Transcript", methods=['GET', 'POST'])
def Transcript():

    # take organism, gene name and enst id from the url
    organism = request.args.get('organism')
    genename = request.args.get('gene')
    enstid = request.args.get('enst')


    if organism == 'Human':
        # Query interactiondisruptionindominanttranscripts from mysql table
        cur2 = mysql.get_db().cursor()
        cur_disease = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId,
                    DomCancerTrans, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_int_anno_Human
                    LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_int_anno_Human.DomCancerTrans = mart_export.`Ensembl Transcript ID`
                    WHERE interactionDisruptionInDominantTranscripts_int_anno_Human.DomCancerTrans = %s '''

        sql_disease = ''' SELECT Tissue, ENSG, NumberOfGtexMDT, GeneName1, GeneName2, TotalNumberOfStringInt, cMDTnumberOfUniqMissedInt, Pfam1, Domain1, Pfam2, Domain2, RNAseqAliquotID,
                    cMDT, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_Human_Disease
                    LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_Human_Disease.cMDT = mart_export.`Ensembl Transcript ID`
                    WHERE interactionDisruptionInDominantTranscripts_Human_Disease.cMDT = %s '''

        adr = (enstid,)
        cur2.execute(sql, adr)
        cur_disease.execute(sql_disease, adr)

        isonet_tuple = cur2.fetchall()
        isonet_tuple_disease = cur_disease.fetchall()

        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2',
                                                'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans',
                                                'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans',
                                                'StringDensityRank1','Region1', 'Transcript_Name'])

        df_disease = pd.DataFrame(isonet_tuple_disease, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2',
                                                'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans',
                                                'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans',
                                                'StringDensityRank1','Region1', 'Transcript_Name'])


        df = pd.concat([df, df_disease])
        df = df.rename(columns={'ENSG': 'ENSGid', 'NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions','TotalNumberOfStringInt': 'NumberOfStringInt'})
        df['Domain1'].replace({"-":"None"},inplace=True)
        df['Domain2'].replace({"-":"None"},inplace=True)
        df['MissedInteractions'].replace({-1:0},inplace=True)

        #result = df.to_dict(orient='records')


        transcript_name = df[df.DomCancerTrans == enstid].iloc[0,15]

        if genename == "":
            genename = df[df.DomCancerTrans == enstid].iloc[0,3]
        elif genename == None:
            genename = df[df.DomCancerTrans == enstid].iloc[0,3]

        enst_list = list(df[df.GeneName1 == genename]['DomCancerTrans'].unique())

        if enstid in enst_list:
            # Query ClinVar genes from database
            cur4 = mysql.get_db().cursor()
            cur4.execute( '''SELECT GeneID, GeneSymbol FROM ClinVar_annotations ''' )
            clinvar_tuple = cur4.fetchall()
            df_clinvar = pd.DataFrame(clinvar_tuple, columns=['GeneID', 'GeneName'])
            df_clinvar_list = df_clinvar['GeneName'].tolist()
            clinvar = df_clinvar.drop_duplicates()
            clinvar_dict = clinvar.to_dict(orient='records')


            df = df.rename(columns={'Tissue':'CancerType2'})

            df = df.drop_duplicates()
            data_dict = df.to_dict(orient='records')

            #make a table for some statistics
            df['Sum'] = df['MissedInteractions']

            # calculate the interaction lost
            statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'Sum', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans']].drop_duplicates()
            statistics_table_dict = statistic_table.to_dict(orient='records')

            string_score = statistic_table.iloc[0,6]
            string_score = int(statistic_table.iloc[0,6]*100)

            TotalInt = int(list(df['NumberOfStringInt'].unique())[0])
            LostInt = int(list(df['Sum'].unique())[0])

            ### draw pie chart to show which cancers/diseases have this transcript as a MDT
            data = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]],
                                subplot_titles=("%Interaction Lost", "Disease Types*"))


            data.add_trace(go.Pie(labels=df.drop_duplicates(subset=['CancerSampleId', 'CancerType2']).CancerType2), 1, 2)

            data.update_traces(selector=dict(type='pie'), textinfo='label+text',hoverinfo='label+percent',
                                      marker=dict(line=dict(color='#000000', width=4)))


            if statistic_table.iloc[0,2] == 0:

                data.add_trace(go.Pie(values=[100, 0], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,1)
                data.update_traces(selector=dict(type='pie'),textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))

            else:
                data.add_trace(go.Pie(values=[(statistic_table.iloc[0,2]-statistic_table.iloc[0,3])*100/(statistic_table.iloc[0,2]),statistic_table.iloc[0,3]*100/(statistic_table.iloc[0,2])], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,1)
                data.update_traces(selector=dict(type='pie'),textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))



            data.update_layout(showlegend=False, title_font_size=18)
            graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)


            ## Take missed interactions for the interested ENST id from the missed_interactions table (isoform interaction network table)
            cur3 = mysql.get_db().cursor()
            sql3 = ''' SELECT ENSG, ENSP, ENST, MissInts FROM interactionsinisoforms_900 WHERE ENST = %s '''
            adr3 = (enstid,)
            cur3.execute(sql3, adr3)
            missed_interactions_tuple = cur3.fetchall()
            missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'MissInts'])
            Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
            Isoform_Int_Network = pd.concat([missed_interactions, Isoform_Int_Network_splitted], axis=1)

            ENSP_list = list()
            ensp_frame = list()

            ## take existing interactions for the interested ENST id from the missed_interactions table (isoform interaction network table)
            cur_ext_int = mysql.get_db().cursor()
            sql_ext_int = ''' SELECT ENSG, ENSP, ENST, ExistInts FROM interactionsinisoforms_900 WHERE ENST = %s '''
            adr_ext_int = (enstid,)
            cur_ext_int.execute(sql_ext_int, adr_ext_int)
            exist_interactions_tuple = cur_ext_int.fetchall()
            exist_interactions = pd.DataFrame(exist_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'ExistInts'])
            Isoform_Int_Network_splitted_exists = pd.DataFrame(exist_interactions.ExistInts.str.split(':').tolist())

            for eachcolumn in range(3, len(Isoform_Int_Network.iloc[0,:])):
                Isoform_Int_Network.iloc[0,eachcolumn] = str(Isoform_Int_Network.iloc[0,eachcolumn])

                if "ENSP" in Isoform_Int_Network.iloc[0,eachcolumn]:
                    ENSP_list.append(Isoform_Int_Network.iloc[0,eachcolumn])

                else:
                    continue

            ensp_frame = ENSP_list

            ## for the existing ints
            ensp_frame_exists = list()
            ENSP_list_exists = list()

            for eachcolumn in range(3, len(Isoform_Int_Network_splitted_exists.iloc[0,:])):
                Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])

                if "ENSP" in Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn]:
                    ENSP_list_exists.append(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])
                else:
                    continue

            ensp_frame_exists = ENSP_list_exists

            partner_genenames = []

            try:
                placeholders = ','.join(['%s'] * len(ensp_frame))
                cur_ensp = mysql.get_db().cursor()
                cur_ensp.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame))
                ensp_tuple = cur_ensp.fetchall()
                ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
                partner_genenames = list(ensp_genename['GeneName'])
            except:
                pass

            partner_genenames_exists = []

            try:
                placeholders = ','.join(['%s'] * len(ensp_frame_exists))
                cur_ensp_exists = mysql.get_db().cursor()
                cur_ensp_exists.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame_exists))
                ensp_tuple_exists = cur_ensp_exists.fetchall()
                ensp_genename_exists = pd.DataFrame(ensp_tuple_exists, columns=['ENSPid', 'GeneName'])
                partner_genenames_exists = list(ensp_genename_exists['GeneName'])
            except:
                pass


            clinvar_partners = [eachpartner for eachpartner in partner_genenames if eachpartner in df_clinvar_list]
            if genename in df_clinvar_list:
                clinvar_partners.append(genename)
            clinvar_partners_df = pd.DataFrame({'GeneName':clinvar_partners})
            clinvar_dict2 = clinvar[clinvar.GeneName.isin(clinvar_partners_df.GeneName)].to_dict(orient='records')


    if organism == 'Mouse':

        # Query interactiondisruptionindominanttranscripts from mysql table
        cur2 = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, NumberOfGTExMDT, GeneName1, GeneName2, TotalNumberOfSTRINGint, cMDTnumberOfUniqMissedInt, NumberOfCommonMissedInt, Pfam1, Domain1, Pfam2, Domain2, RNAseqAliquotID,
                    cMDT, StringDensityRank1, Region1, mart_export_mouse.Transcript_Name FROM interactionDisruptionInDominantTranscripts_Mouse
                    LEFT JOIN mart_export_mouse ON interactionDisruptionInDominantTranscripts_Mouse.cMDT = mart_export_mouse.ENSMUST
                    WHERE interactionDisruptionInDominantTranscripts_Mouse.cMDT = %s '''
        adr = (enstid,)
        cur2.execute(sql, adr)
        isonet_tuple = cur2.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2',
                                                'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'NumberOfCommonMissedInt',
                                                'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans',
                                                'StringDensityRank1','Region1', 'Transcript_Name'])

        df = df.rename(columns={'Tissue': 'CancerType2', 'ENSG': 'ENSGid', 'NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions','TotalNumberOfStringInt': 'NumberOfStringInt'})
        df['Domain1'].replace({"-":"None"},inplace=True)
        df['Domain2'].replace({"-":"None"},inplace=True)
        df['MissedInteractions'].replace({-1:0},inplace=True)
        df['NumberOfCommonMissedInt'].replace({-1:0},inplace=True)
        df['Sum'] = df['MissedInteractions'] + df['NumberOfCommonMissedInt']

        transcript_name = df[df.DomCancerTrans == enstid].iloc[0,16]

        if genename == "":
            genename = df[df.DomCancerTrans == enstid].iloc[0,3]
        elif genename == None:
            genename = df[df.DomCancerTrans == enstid].iloc[0,3]

        enst_list = list(df[df.GeneName1 == genename]['DomCancerTrans'].unique())

        if enstid in enst_list:
            # Query ClinVar genes from database
            cur4 = mysql.get_db().cursor()
            cur4.execute( '''SELECT GeneID, GeneSymbol FROM ClinVar_annotations ''' )
            clinvar_tuple = cur4.fetchall()
            df_clinvar = pd.DataFrame(clinvar_tuple, columns=['GeneID', 'GeneName'])
            df_clinvar_list = df_clinvar['GeneName'].tolist()
            clinvar = df_clinvar.drop_duplicates()
            clinvar_dict = clinvar.to_dict(orient='records')

            data_dict = df.to_dict(orient='records')

            #calculate the percentage of interaction lost
            statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'Sum', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans']].drop_duplicates()
            statistics_table_dict = statistic_table.to_dict(orient='records')

            TotalInt = int(list(df['NumberOfStringInt'].unique())[0])
            LostInt = int(list(df['Sum'].unique())[0])

            string_score = statistic_table.iloc[0,6]
            string_score = int(statistic_table.iloc[0,6]*100)

            ### draw pie chart to show which cancers/diseases have this transcript as a MDT
            data = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]],
                                subplot_titles=("% Interaction Lost", "Disease Types"))


            data.add_trace(go.Pie(labels=df.drop_duplicates(subset=['CancerSampleId', 'CancerType2']).CancerType2), 1, 2)

            data.update_traces(selector=dict(type='pie'), textinfo='label+text',hoverinfo='label+percent',
                                      marker=dict(line=dict(color='#000000', width=4)))


            if statistic_table.iloc[0,2] == 0:

                data.add_trace(go.Pie(values=[100, 0], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,1)
                data.update_traces(selector=dict(type='pie'),textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))

            else:
                data.add_trace(go.Pie(values=[(statistic_table.iloc[0,2]-statistic_table.iloc[0,3])*100/(statistic_table.iloc[0,2]),statistic_table.iloc[0,3]*100/(statistic_table.iloc[0,2])], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,1)
                data.update_traces(selector=dict(type='pie'),textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))



            data.update_layout(showlegend=False, title_font_size=18)
            graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)


            ## Take missed interactions for the interested ENST id from the missed_interactions table (isoform specific interaction table)
            cur3 = mysql.get_db().cursor()
            sql3 = ''' SELECT ENSG, ENSP, ENST, MissInts FROM interactionsinisoforms_mouse WHERE ENST = %s '''
            adr3 = (enstid,)
            cur3.execute(sql3, adr3)
            missed_interactions_tuple = cur3.fetchall()
            missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'MissInts'])
            Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
            Isoform_Int_Network = pd.concat([missed_interactions, Isoform_Int_Network_splitted], axis=1)

            ENSP_list = list()
            ensp_frame = list()

            ## take existing interactions for the interested ENST id from the missed_interactions table (isoform specific interaction table)
            cur_ext_int = mysql.get_db().cursor()
            sql_ext_int = ''' SELECT ENSG, ENSP, ENST, ExistInts FROM interactionsinisoforms_mouse WHERE ENST = %s '''
            adr_ext_int = (enstid,)
            cur_ext_int.execute(sql_ext_int, adr_ext_int)
            exist_interactions_tuple = cur_ext_int.fetchall()
            exist_interactions = pd.DataFrame(exist_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'ExistInts'])
            Isoform_Int_Network_splitted_exists = pd.DataFrame(exist_interactions.ExistInts.str.split(':').tolist())

            for eachcolumn in range(3, len(Isoform_Int_Network.iloc[0,:])):
                Isoform_Int_Network.iloc[0,eachcolumn] = str(Isoform_Int_Network.iloc[0,eachcolumn])

                if "ENSMUSP" in Isoform_Int_Network.iloc[0,eachcolumn]:
                    ENSP_list.append(Isoform_Int_Network.iloc[0,eachcolumn])

                else:
                    continue

            ensp_frame = ENSP_list

            ## for the existing ints
            ensp_frame_exists = list()
            ENSP_list_exists = list()

            for eachcolumn in range(3, len(Isoform_Int_Network_splitted_exists.iloc[0,:])):
                Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])

                if "ENSMUSP" in Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn]:
                    ENSP_list_exists.append(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])
                else:
                    continue

            ensp_frame_exists = ENSP_list_exists

            partner_genenames = []

            try:
                placeholders = ','.join(['%s'] * len(ensp_frame))
                cur_ensp = mysql.get_db().cursor()
                cur_ensp.execute('''SELECT ENSPid, GeneName FROM ensg_enst_ensp_des_mouse WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame))
                ensp_tuple = cur_ensp.fetchall()
                ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
                partner_genenames = list(ensp_genename['GeneName'])
            except:
                pass

            partner_genenames_exists = []

            try:
                placeholders = ','.join(['%s'] * len(ensp_frame_exists))
                cur_ensp_exists = mysql.get_db().cursor()
                cur_ensp_exists.execute('''SELECT ENSPid, GeneName FROM ensg_enst_ensp_des_mouse WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame_exists))
                ensp_tuple_exists = cur_ensp_exists.fetchall()
                ensp_genename_exists = pd.DataFrame(ensp_tuple_exists, columns=['ENSPid', 'GeneName'])
                partner_genenames_exists = list(ensp_genename_exists['GeneName'])
            except:
                pass


            clinvar_partners = [eachpartner for eachpartner in partner_genenames if eachpartner in df_clinvar_list]
            if genename in df_clinvar_list:
                clinvar_partners.append(genename)
            clinvar_partners_df = pd.DataFrame({'GeneName':clinvar_partners})
            clinvar_dict2 = clinvar[clinvar.GeneName.isin(clinvar_partners_df.GeneName)].to_dict(orient='records')

    return render_template('network.html', ensp_frame_exists=ensp_frame_exists, string_score=string_score, df_clinvar_list=df_clinvar_list, partner_genenames_exists=partner_genenames_exists, transcript_name=transcript_name, genename=genename, organism = organism, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, graphJSON2=graphJSON2, clinvar_dict=clinvar_dict, clinvar_dict2=clinvar_dict2, TotalInt = TotalInt, LostInt=LostInt)



### this page is for diseases excluding cancer
@app.route("/Disease", methods=['GET', 'POST'])
def Disease():

    SampleDiseaseType = request.args.get('disease')
    # mouse diseases
    if (SampleDiseaseType == 'OIH_NAc' or SampleDiseaseType == 'OIH_TG'):
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, GeneName1, RNAseqAliquotID, cMDT, cMDTrelExpression, mart_export_mouse.Transcript_Name FROM interactionDisruptionInDominantTranscripts_Mouse
        LEFT JOIN mart_export_mouse ON interactionDisruptionInDominantTranscripts_Mouse.cMDT = mart_export_mouse.ENSMUST
        WHERE interactionDisruptionInDominantTranscripts_Mouse.Tissue = %s '''
        adr = (SampleDiseaseType, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'GeneName1', 'RNAseqAliquotID', 'cMDT',  'dMDTenrichment', 'Transcript_Name'])
        df_iso2 = df.drop_duplicates()
        df_iso2 = df_iso2[['Tissue', 'GeneName1', 'RNAseqAliquotID', 'cMDT', 'Transcript_Name', 'dMDTenrichment']]
        temp_dict = df_iso2.to_dict(orient='records')

        ## second table in the page representing second graph - MDT counts per sample
        cur2 = mysql.get_db().cursor()
        sql2 = '''SELECT Tissue, SampleID, Count FROM sample_mdt_counts WHERE Tissue = %s '''
        adr2 = (SampleDiseaseType, )
        cur2.execute(sql2, adr2)
        dataset_muts = cur2.fetchall()
        muts_df = pd.DataFrame(dataset_muts, columns= ['Tissue', 'SampleID', 'Count' ])
        dataset = muts_df.to_dict(orient='records')
    else:
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, GeneName1, RNAseqAliquotID, cMDT, cMDTrelExpression, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_Human_Disease
        LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_Human_Disease.cMDT = mart_export.`Ensembl Transcript ID`
        WHERE interactionDisruptionInDominantTranscripts_Human_Disease.Tissue = %s '''
        adr = (SampleDiseaseType, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'GeneName1', 'RNAseqAliquotID', 'cMDT', 'dMDTenrichment', 'Transcript_Name'])
        df_iso2 = df.drop_duplicates()
        df_iso2 = df_iso2[['Tissue', 'GeneName1', 'RNAseqAliquotID', 'cMDT', 'Transcript_Name', 'dMDTenrichment']]
        temp_dict = df_iso2.to_dict(orient='records')

        ## second table in the page representing second graph - MDT counts per sample
        cur2 = mysql.get_db().cursor()
        sql2 = '''SELECT Tissue, SampleID, Count FROM sample_mdt_counts WHERE Tissue = %s '''
        adr2 = (SampleDiseaseType, )
        cur2.execute(sql2, adr2)
        dataset_muts = cur2.fetchall()
        muts_df = pd.DataFrame(dataset_muts, columns= ['Tissue', 'SampleID', 'Count' ])
        dataset = muts_df.to_dict(orient='records')

    return render_template("Disease_Based.html", SampleDiseaseType = SampleDiseaseType, data = temp_dict, data2 = dataset)

def DiseaseSpecific(SampleDiseaseType):

    SampleDiseaseType = request.args.get('disease')

    if (SampleDiseaseType == 'OIH_NAc'or SampleDiseaseType == 'OIH_TG') :
        cur = mysql.get_db().cursor()
        sql = ''' SELECT cMDT, Frequency, Tissue, mart_export_mouse.Transcript_Name FROM dMDT_Frequency LEFT JOIN mart_export_mouse ON dMDT_Frequency.cMDT = mart_export_mouse.ENSMUST
                WHERE dMDT_Frequency.Tissue = %s ORDER BY Frequency DESC LIMIT 10 '''
        adr = (SampleDiseaseType, )
        cur.execute(sql,adr)
        TableS4_tuple = cur.fetchall()
        result = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','Tissue', 'Transcript_Name'])

    else:
        cur = mysql.get_db().cursor()
        sql = ''' SELECT cMDT, Frequency, Tissue, mart_export.`Associated Transcript Name` FROM dMDT_Frequency LEFT JOIN mart_export ON dMDT_Frequency.cMDT = mart_export.`Ensembl Transcript ID`
            WHERE dMDT_Frequency.Tissue = %s ORDER BY Frequency DESC LIMIT 10 '''
        adr = (SampleDiseaseType, )
        cur.execute(sql,adr)
        TableS4_tuple = cur.fetchall()
        result = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','Tissue', 'Transcript_Name'])


    result['Frequency'] = result['Frequency'].replace(np.nan, 0)
    result['Frequency'] = pd.to_numeric(result["Frequency"])

    cur2 = mysql.get_db().cursor()
    sql2 = ''' SELECT Tissue, SampleID, Count  FROM sample_mdt_counts WHERE Tissue = %s '''
    adr2 = (SampleDiseaseType, )
    cur2.execute(sql2, adr2)
    cMDT_mdts_muts = cur2.fetchall()
    cMDT_dist = pd.DataFrame(cMDT_mdts_muts, columns=['Tissue', 'SampleID','Count'])
    cMDT_dist['Count'] = pd.to_numeric(cMDT_dist["Count"])



    # extract dMDT TPMs values from pcawg data - Top 10 dMDT
    enst_list = result['cMDT'].values.tolist()
    tissue = request.args.get('disease')
    pcawg = '_pcawg'

    enst_tpms = mysql.get_db().cursor()
    my_tissue2 = tissue + pcawg
    result_enst_tpms = pd.DataFrame()
    for eachdMDT in enst_list:
        enst_tpms_sql = ''' SELECT * FROM `''' + my_tissue2 + '''` WHERE Feature REGEXP %s '''
        enst_tpms_adr = (eachdMDT,)
        enst_tpms.execute(enst_tpms_sql, enst_tpms_adr)
        enst_tpms_tuple = enst_tpms.fetchall()
        result_enst_tpms_df = pd.DataFrame(enst_tpms_tuple)
        result_enst_tpms = result_enst_tpms.append(result_enst_tpms_df, ignore_index = True)

    medians_each_row_list = []
    result_enst_tpms_TPMs = result_enst_tpms.iloc[:,1:]

    medians_each_row = result_enst_tpms_TPMs.median(axis=1)

    medians_each_row_list = medians_each_row.values.tolist()
    medians_rounds = [round(number, 2) for number in medians_each_row_list]
    medians_rounds_text = ["Median TPM: " + str(item) for item in medians_rounds]

    plot = make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3], subplot_titles=("Frequency of top 10 dMDT with median TPM values", "Number of dMDT"))




    trace1 = go.Bar(
                #name = '% of Transcripts across ' + SampleDiseaseType,
                marker_color = 'indianred',
                x=result.iloc[:,3],
                y=result.iloc[:,1]*100,
                text = medians_rounds,
                hovertext = medians_rounds_text,
                textposition='inside'
            )

    trace2 = go.Box(y=cMDT_dist.Count, boxpoints='outliers',jitter=0.3,
                    marker_color = '#00CC96', name=' ')


    plot.append_trace(trace1, row=1, col=1)
    plot.update_yaxes(title_text="Frequency of Transcripts in Samples (%)", row=1, col=1)
    plot.update_yaxes(title_text="Distribiton of dMDT Across Samples",row=1, col=2)
    plot.update_xaxes(title_text=' ', row=1, col=2)
    plot.append_trace(trace2, row=1, col=2)
    plot.update_layout(showlegend=False)

    graphJSON = json.dumps(plot, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/DiseaseSpecific', methods=['GET', 'POST'])
def change_features_disease():
    SampleDiseaseType = request.args.get('disease')
    graphJSON = DiseaseSpecific(SampleDiseaseType)
    return graphJSON


@app.route("/Gene", methods=['GET', 'POST'])
def Gene():

    # Take organisim and genename from url
    genename = request.args.get('gene')
    organism = request.args.get('organism')

    if organism == 'Human':
        # Query ClinVar genes from database
        cur4 = mysql.get_db().cursor()
        cur4.execute( '''SELECT GeneID, GeneSymbol FROM ClinVar_annotations ''' )
        clinvar_tuple = cur4.fetchall()
        df_clinvar = pd.DataFrame(clinvar_tuple, columns=['GeneID', 'GeneName'])
        df_clinvar_list = df_clinvar['GeneName'].tolist()
        clinvar = df_clinvar.drop_duplicates()
        clinvar_dict = clinvar.to_dict(orient='records')

        cur5 = mysql.get_db().cursor()
        cur_disease = mysql.get_db().cursor()


        sql5 = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_int_anno_Human
        LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_int_anno_Human.DomCancerTrans = mart_export.`Ensembl Transcript ID`
        WHERE interactionDisruptionInDominantTranscripts_int_anno_Human.GeneName1 = %s '''

        sql_disease = ''' SELECT Tissue, ENSG, NumberOfGTExMDT, GeneName1, GeneName2, TotalNumberOfStringInt, cMDTnumberOfUniqMissedInt, Pfam1, Domain1, Pfam2, Domain2, RNAseqAliquotID, cMDT, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactionDisruptionInDominantTranscripts_Human_Disease
        LEFT JOIN mart_export ON interactionDisruptionInDominantTranscripts_Human_Disease.cMDT = mart_export.`Ensembl Transcript ID`
        WHERE interactionDisruptionInDominantTranscripts_Human_Disease.GeneName1 = %s '''

        adr5 = (genename,)

        cur5.execute(sql5, adr5)
        isonet_tuple = cur5.fetchall()

        cur_disease.execute(sql_disease, adr5)
        isonet_tuple_disease = cur_disease.fetchall()

        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'cMDTnumberOfUniqMissedInt', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1', 'Region1', 'Transcript_Name'])
        df =  df.rename(columns={'Tissue': 'CancerType2', 'ENSG': 'ENSGid', 'TotalNumberOfStringInt': 'NumberOfStringInt','cMDTnumberOfUniqMissedInt': 'MissedInteractions'})

        df_disease = pd.DataFrame(isonet_tuple_disease, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'cMDTnumberOfUniqMissedInt', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1', 'Region1', 'Transcript_Name'])
        df_disease =  df_disease.rename(columns={'Tissue': 'CancerType2', 'ENSG': 'ENSGid', 'TotalNumberOfStringInt': 'NumberOfStringInt','cMDTnumberOfUniqMissedInt': 'MissedInteractions'})

        df = pd.concat([df, df_disease])

        #df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
        df = df.drop_duplicates()

        df['Domain1'].replace({"-":"None"},inplace=True)
        df['Domain2'].replace({"-":"None"},inplace=True)
        df['MissedInteractions'].replace({-1:0},inplace=True)

        result = df.to_dict(orient='records')

        statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans', 'Transcript_Name']].drop_duplicates()
        statistics_table_dict = statistic_table.to_dict(orient='records')

        data_dict = df.to_dict(orient='records')

        string_score = statistic_table.iloc[0,6]
        string_score = int(statistic_table.iloc[0,6]*100)
        #string_score = float("{:.2f}".format(string_score))*100
        ### DRAW PIE CHART

        data = make_subplots(rows=1, cols=1, specs=[[{'type':'domain'}]])

        data.add_trace(go.Pie(labels=df.CancerType2, title="Disease Types"), 1, 1)
        data.update_traces(selector=dict(type='pie'), textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(line=dict(color='#000000', width=4)))

        data.update_layout(title=" Gene Name: {} ".format(genename), showlegend=False, title_font_size=24)
        graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

        # take missed interactions from the db
        ensg = df.iloc[0,1]
        cur3 = mysql.get_db().cursor()
        sql3 = ''' SELECT ENSG, ENSP, ENST, ExistInts, MissInts FROM interactionsinisoforms_900 WHERE ENSG = %s LIMIT 1'''
        adr3 = (ensg,)
        cur3.execute(sql3, adr3)
        missed_interactions_tuple = cur3.fetchall()
        missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'ExistsInts', 'MissInts'])
        Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
        Isoform_Int_Network_splitted_exists = pd.DataFrame(missed_interactions.ExistsInts.str.split(':').tolist())

        ENSP_list = list()
        ensp_frame = list()


        for eachcolumn in range(4, len(Isoform_Int_Network_splitted.iloc[0,:])):
            Isoform_Int_Network_splitted.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted.iloc[0,eachcolumn])

            if "ENSP" in Isoform_Int_Network_splitted.iloc[0,eachcolumn]:
                ENSP_list.append(Isoform_Int_Network_splitted.iloc[0,eachcolumn])

            else:
                continue

        ensp_frame = ENSP_list
        ensp_frame_exists = list()
        ENSP_list_exists = list()


        for eachcolumn in range(3, len(Isoform_Int_Network_splitted_exists.iloc[0,:])):
            Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])

            if "ENSP" in Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn]:
                ENSP_list_exists.append(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])
            else:
                continue

        ensp_frame_exists = ENSP_list_exists

        partner_genenames = []

        try:
            placeholders = ','.join(['%s'] * len(ensp_frame))
            cur_ensp = mysql.get_db().cursor()
            cur_ensp.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame))
            ensp_tuple = cur_ensp.fetchall()
            ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
            partner_genenames = list(ensp_genename['GeneName'])
        except:
            pass

        partner_genenames_exists = []

        try:
            placeholders = ','.join(['%s'] * len(ensp_frame_exists))
            cur_ensp_exists = mysql.get_db().cursor()
            cur_ensp_exists.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame_exists))
            ensp_tuple_exists = cur_ensp_exists.fetchall()
            ensp_genename_exists = pd.DataFrame(ensp_tuple_exists, columns=['ENSPid', 'GeneName'])
            partner_genenames_exists = list(ensp_genename_exists['GeneName'])
        except:
            pass

    if organism == 'Mouse':

        cur4 = mysql.get_db().cursor()
        cur4.execute( '''SELECT GeneID, GeneSymbol FROM ClinVar_annotations ''' )
        clinvar_tuple = cur4.fetchall()
        df_clinvar = pd.DataFrame(clinvar_tuple, columns=['GeneID', 'GeneName'])
        df_clinvar_list = df_clinvar['GeneName'].tolist()
        clinvar = df_clinvar.drop_duplicates()
        clinvar_dict = clinvar.to_dict(orient='records')

        cur5 = mysql.get_db().cursor()
        sql5 = ''' SELECT Tissue, ENSG, NumberOfGTExMDT, GeneName1, GeneName2, TotalNumberOfStringInt, cMDTnumberOfUniqMissedInt, Pfam1, Domain1, Pfam2, Domain2, RNAseqAliquotID, cMDT, StringDensityRank1, Region1, mart_export_mouse.Transcript_Name FROM interactionDisruptionInDominantTranscripts_Mouse
        LEFT JOIN mart_export_mouse ON interactionDisruptionInDominantTranscripts_Mouse.cMDT = mart_export_mouse.ENSMUST
        WHERE interactionDisruptionInDominantTranscripts_Mouse.GeneName1 = %s '''
        adr5 = (genename,)
        cur5.execute(sql5, adr5)
        isonet_tuple = cur5.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1', 'Region1', 'Transcript_Name'])
        df =  df.rename(columns={'Tissue': 'CancerType2','ENSG': 'ENSGid', 'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions'})
        df['Domain1'].replace({"-":"None"},inplace=True)
        df['Domain2'].replace({"-":"None"},inplace=True)
        df['MissedInteractions'].replace({-1:0},inplace=True)

        result = df.to_dict(orient='records')

        statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans', 'Transcript_Name']].drop_duplicates()
        statistics_table_dict = statistic_table.to_dict(orient='records')

        data_dict = df.to_dict(orient='records')

        string_score = statistic_table.iloc[0,6]
        string_score = int(statistic_table.iloc[0,6]*100)
        #string_score = float("{:.2f}".format(string_score))*100
        ### DRAW PIE CHART

        data = make_subplots(rows=1, cols=1, specs=[[{'type':'domain'}]])

        data.add_trace(go.Pie(labels=df.CancerType2, title="Cancer Types"), 1, 1)
        data.update_traces(selector=dict(type='pie'), textinfo='label+text',hoverinfo='label+percent',
                                  marker=dict(line=dict(color='#000000', width=4)))

        data.update_layout(title=" Gene Name: {} ".format(genename), showlegend=False, title_font_size=24)
        graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

        # take missed interactions from the db
        ensg = df.iloc[0,1]
        cur3 = mysql.get_db().cursor()
        sql3 = ''' SELECT ENSG, ENSP, ENST, ExistInts, MissInts FROM interactionsinisoforms_mouse WHERE ENSG = %s LIMIT 1'''
        adr3 = (ensg,)
        cur3.execute(sql3, adr3)
        missed_interactions_tuple = cur3.fetchall()
        missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'ExistsInts', 'MissInts'])
        Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
        #Isoform_Int_Network = pd.concat([missed_interactions, Isoform_Int_Network_splitted], axis=1)

        Isoform_Int_Network_splitted_exists = pd.DataFrame(missed_interactions.ExistsInts.str.split(':').tolist())

        ENSP_list = list()
        ensp_frame = list()


        for eachcolumn in range(4, len(Isoform_Int_Network_splitted.iloc[0,:])):
            Isoform_Int_Network_splitted.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted.iloc[0,eachcolumn])

            if "ENSP" in Isoform_Int_Network_splitted.iloc[0,eachcolumn]:
                ENSP_list.append(Isoform_Int_Network_splitted.iloc[0,eachcolumn])

            else:
                continue

        ensp_frame = ENSP_list
        ensp_frame_exists = list()
        ENSP_list_exists = list()


        for eachcolumn in range(3, len(Isoform_Int_Network_splitted_exists.iloc[0,:])):
            Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn] = str(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])

            if "ENSMUSP" in Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn]:
                ENSP_list_exists.append(Isoform_Int_Network_splitted_exists.iloc[0,eachcolumn])
            else:
                continue

        ensp_frame_exists = ENSP_list_exists

        partner_genenames = []

        try:
            placeholders = ','.join(['%s'] * len(ensp_frame))
            cur_ensp = mysql.get_db().cursor()
            cur_ensp.execute('''SELECT ENSPid, GeneName FROM ensg_enst_ensp_des_mouse WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame))
            ensp_tuple = cur_ensp.fetchall()
            ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
            partner_genenames = list(ensp_genename['GeneName'])
        except:
            pass

        partner_genenames_exists = []

        try:
            placeholders = ','.join(['%s'] * len(ensp_frame_exists))
            cur_ensp_exists = mysql.get_db().cursor()
            cur_ensp_exists.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des_mouse WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame_exists))
            ensp_tuple_exists = cur_ensp_exists.fetchall()
            ensp_genename_exists = pd.DataFrame(ensp_tuple_exists, columns=['ENSPid', 'GeneName'])
            partner_genenames_exists = list(ensp_genename_exists['GeneName'])
        except:
            pass

    return render_template('GeneBased.html', genename=genename, organism = organism, partner_genenames_exists=partner_genenames_exists, df_clinvar_list=df_clinvar_list, partner_genenames=partner_genenames, string_score=string_score, data=data_dict, data_statistics = statistics_table_dict, result=result, graphJSON2=graphJSON2, clinvar_dict=clinvar_dict)


@app.route("/Sample", methods=['GET', 'POST'])
def Sample():

    # get sample id, gene name and tissue type from url
    CancerSampleId = request.args.get('sampleid')
    genename = request.args.get('gene')
    tissue = request.args.get('tissue')

    if (tissue == 'OIH_NAc' or tissue == 'OIH_TG'):
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, RNAseqAliquotID, cMDT, GTExMDT FROM interactionDisruptionInDominantTranscripts_Mouse WHERE RNAseqAliquotID = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])


    elif (tissue == "Alzheimer's disease" or tissue == "Parkinson's disease" or tissue == "Behcet's disease"):
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, RNAseqAliquotID, cMDT, GTExMDT FROM interactionDisruptionInDominantTranscripts_Human_Disease WHERE RNAseqAliquotID = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])

    else:
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactionDisruptionInDominantTranscripts_int_anno_Human WHERE CancerSampleId = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])

    df =  df.rename(columns={'ENSG': 'ENSGid'})
    df = df.drop_duplicates()

    dff = df[df.GeneName1 == genename]
    cancer_trans_id = list(dff['DomCancerTrans'])[0] ## DomCancerTrans id

    normal_trans_ids = dff['GTExMDIs'].str.split(";", expand = True) # n=1
    normal_trans_id_list = []

    #list of normal trans id
    for i in range(len(normal_trans_ids.iloc[0,:])):
        for j in range(len(normal_trans_ids.iloc[:,0])):
            if "ENST" in str(normal_trans_ids.iloc[j,i]):
                normal_trans_id_list.append(normal_trans_ids.iloc[j,i].split(":", 1)[0])

            if "ENSMUST" in str(normal_trans_ids.iloc[j,i]):
                normal_trans_id_list.append(normal_trans_ids.iloc[j,i].split(":", 1)[0])
            else:
                continue

    normal_trans_id_dict = dict()
    for index,value in enumerate(normal_trans_id_list):
      normal_trans_id_dict[index] = value

    return render_template("Sample_Based.html", CancerSampleId=CancerSampleId, tissue=tissue, genename=genename, normal_trans_id_dict=normal_trans_id_dict, normal_trans_id_list=normal_trans_id_list, cancer_trans_id=cancer_trans_id)

def update_fig(CancerSampleId, genename, tissue):
    tissue = request.args.get('tissuetype')
    CancerSampleId = request.args.get('CanSampleId')

    if (tissue == 'OIH_NAc' or tissue == 'OIH_TG'):
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, RNAseqAliquotID, cMDT, GTExMDT FROM interactionDisruptionInDominantTranscripts_Mouse WHERE RNAseqAliquotID = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
        tissuetype = tissue
    elif (tissue == "Alzheimer's disease" or tissue == "Parkinson's disease" or tissue == "Behcet's disease"):
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, RNAseqAliquotID, cMDT, GTExMDT FROM interactionDisruptionInDominantTranscripts_Human_Disease WHERE RNAseqAliquotID = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
        tissuetype = tissue
    else:
        cur = mysql.get_db().cursor()
        sql = ''' SELECT Tissue, ENSG, GeneName1, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactionDisruptionInDominantTranscripts_int_anno_Human WHERE CancerSampleId = %s '''
        adr = (CancerSampleId, )
        cur.execute(sql, adr)
        isonet_tuple = cur.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
        tissuetype = tissue.split('.')[1].replace('-','_')

    df =  df.rename(columns={'ENSG': 'ENSGid'})
    dff = df[df.GeneName1 == genename]
    cancer_trans_id = dff.iloc[0,4] ## DomCancerTrans id
    normal_trans_ids = dff['GTExMDIs'].str.split(";", expand = True) # n=1
    normal_trans_id_list = []

     #list of normal trans id

    for i in range(len(normal_trans_ids.iloc[0,:])):
        for j in range(len(normal_trans_ids.iloc[:,0])):
            if "ENST" in str(normal_trans_ids.iloc[j,i]):
                normal_trans_id_list.append(normal_trans_ids.iloc[j,i].split(":", 1)[0])
            elif "ENSMUST" in str(normal_trans_ids.iloc[j,i]):
                normal_trans_id_list.append(normal_trans_ids.iloc[j,i].split(":", 1)[0])
            else:
                continue


    # extract normal transcripts expressions from gtex data - it can be more than 1 transcript
    gtex = '_gtex'
    my_tissue = tissuetype+gtex
    gtex_cur_normal = mysql.get_db().cursor()
    df_gtex_normal = pd.DataFrame()
    for eachtranscript in normal_trans_id_list:

        gtex_sql_normal = '''SELECT * FROM `''' + my_tissue + '''`WHERE Feature REGEXP %s'''
        gtex_adr_normal = (eachtranscript,)
        gtex_cur_normal.execute(gtex_sql_normal, gtex_adr_normal)
        df_gtex_tuple = gtex_cur_normal.fetchall()
        df_gtex_normal_df = pd.DataFrame(df_gtex_tuple)
        df_gtex_normal = df_gtex_normal.append(df_gtex_normal_df, ignore_index = True)

    # extract cancer transcript expression from gtex data - only 1 transcript

    gtex_cur_cancer = mysql.get_db().cursor()
    gtex_sql_cancer = ''' SELECT * FROM `''' + my_tissue + '''` WHERE Feature REGEXP %s '''
    gtex_adr_cancer = (cancer_trans_id,)
    gtex_cur_cancer.execute(gtex_sql_cancer, gtex_adr_cancer)
    df_gtex_tuple_2 = gtex_cur_cancer.fetchall()
    df_gtex_cancer = pd.DataFrame(df_gtex_tuple_2)

    pcawg = '_pcawg'
    my_tissue2 = tissuetype + pcawg


    #extract normal transcripts expressions from pcawg data - it can be more than 1 transcript

    pcawg_cur_normal = mysql.get_db().cursor()
    df_pcawg_normal = pd.DataFrame()
    for eachtranscript in normal_trans_id_list:
        pcawg_sql_normal = ''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM `''' + my_tissue2 + '''` WHERE Feature REGEXP %s '''
        pcawg_adr_normal = (eachtranscript,)
        pcawg_cur_normal.execute(pcawg_sql_normal, pcawg_adr_normal)
        df_pcawg_tuple = pcawg_cur_normal.fetchall()
        df_pcawg_normal_df = pd.DataFrame(df_pcawg_tuple)
        df_pcawg_normal = df_pcawg_normal.append(df_pcawg_normal_df, ignore_index = True)

    # extract cancer transcript expression from pcawg data - only 1 transcript
    pcawg_cur_cancer = mysql.get_db().cursor()
    pcawg_sql_cancer = ''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM `''' + my_tissue2 + '''` WHERE Feature REGEXP %s '''
    pcawg_adr_cancer = (cancer_trans_id,)
    pcawg_cur_cancer.execute(pcawg_sql_cancer, pcawg_adr_cancer)
    df_pcawg_tuple_2 = pcawg_cur_cancer.fetchall()
    df_pcawg_cancer = pd.DataFrame(df_pcawg_tuple_2)

    data = go.Figure()

    ## Transcript Counts of Transcripts in gtex data
    std_gtex_cancer = statistics.pstdev(df_gtex_cancer.iloc[0,1:len(df_gtex_cancer.columns)])
    count_data_cancer_trans_median = np.median(df_gtex_cancer.iloc[0,1:len(df_gtex_cancer.columns)])
    std_gtex_normal= []
    count_data_normal_trans_median = []
    count_df = pd.DataFrame()

    #for each_normal in list(df_gtex_normal.iloc[:,0]):
    for each_normal in set(list(df_gtex_normal.iloc[:,0])):
        x = [cancer_trans_id]
        x.append(each_normal)

        count_data_normal_trans_median = np.median(df_gtex_normal[df_gtex_normal.iloc[:,0] == each_normal].iloc[0,1:len(df_gtex_normal.columns)])
        std_gtex_normal = statistics.pstdev(df_gtex_normal[df_gtex_normal.iloc[:,0] == each_normal].iloc[0,1:len(df_gtex_normal.columns)])

        count_df = count_df.append(pd.DataFrame({"ENST": each_normal,"Normal":count_data_normal_trans_median, "Disease":float(df_pcawg_normal[df_pcawg_normal.iloc[:,0] == each_normal].iloc[0,1]), "Std":std_gtex_normal}, index=[each_normal]))

    cancer_exp = pd.DataFrame({"ENST": cancer_trans_id, "Normal":count_data_cancer_trans_median, "Disease":float(df_pcawg_cancer.iloc[:,1]), "Std":std_gtex_cancer}, index=[cancer_trans_id])
    count_df = count_df.append(cancer_exp)

    for col in count_df.iloc[:,1:2].columns:
        data.add_trace(go.Bar(x=count_df.index, y=count_df[col], name = col, error_y=dict(type='data', array=list(count_df['Std'])),marker_color='indianred'))

    for col in count_df.iloc[:,2:3].columns:
        data.add_trace(go.Bar(x=count_df.index, y=count_df[col], name = col, marker_color='lightsalmon'))

    data.update_layout(title=genename, yaxis_title="TPM Count", title_font_size=24)

    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/SampleBased', methods=['GET', 'POST'])
def change_features4():

    CancerSampleId = request.args.get('CanSampleId')
    genename = request.args.get('genename')
    tissue = request.args.get('tissuetype')

    graphJSON = update_fig(CancerSampleId, genename, tissue)

    return graphJSON


@app.errorhandler(500)
def page_not_found(e):

    genename = request.args.get('gene')
    enstid = request.args.get('enst')

    return render_template('500.html', genename=genename, enstid=enstid)

if __name__ == "__main__":
  app.run(debug=True)