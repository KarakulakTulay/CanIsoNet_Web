from flask import Flask, render_template, request
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import json
from flaskext.mysql import MySQL
import functools
import operator
import statistics
import yaml

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

    cur = mysql.get_db().cursor()
    cur.execute('''SELECT * FROM cancertypesinproject''')
    cancer_types = cur.fetchall()
    cancertypes = pd.DataFrame(cancer_types,columns=['PCAWG-Code', 'Cancer Type Name', 'PCAWG_GTEx'])
   # cancertypes_dict = list(cancertypes.iloc[:,2])
    cancertypes_dict = cancertypes.to_dict(orient='records')

    cur2 = mysql.get_db().cursor()
    cur2.execute(''' SELECT * FROM ENST_Genename_ENSG_TranscriptName ''')
    genename_cmdt = cur2.fetchall()
    genenamecmdt = pd.DataFrame(genename_cmdt, columns=['Feature', 'DomCancerTrans','GeneName1_x', 'ENSG', 'Transcript_Name'])
    genenamecmdt_gene_list = list(genenamecmdt.iloc[:,2].unique())
    genenamecmdt_gene_cmdt = list(genenamecmdt.iloc[:,1].unique())
    genenamecmdt_gene_tn = list(genenamecmdt.iloc[:,4].unique())
    temp_dict = genenamecmdt.to_dict(orient='records')

    return render_template('home.html', data=cancertypes_dict, data2=genenamecmdt_gene_list, data3=genenamecmdt_gene_cmdt, temp_dict=temp_dict, genenamecmdt_gene_tn=genenamecmdt_gene_tn)


@app.route("/cspec_cMDTs")
def cspec_cMDTs():

    cur = mysql.get_db().cursor()
    cur.execute(''' SELECT CancerType, cMDT, GeneName, Count, Total, Frequency, mart_export.`Associated Transcript Name` FROM tables9 LEFT JOIN mart_export ON tables9.cMDT = mart_export.`Ensembl Transcript ID` ORDER BY Frequency DESC''')
    TableS9_tuple = cur.fetchall()
    result = pd.DataFrame(TableS9_tuple, columns=['CancerType', 'cMDT', 'GeneName','Count',  'Total', 'Frequency', 'Transcript_Name'])
    result[['Splitted','CancerType2']] = result.CancerType.str.split('.', expand=True)
    result2 = result[['CancerType2', 'CancerType', 'cMDT', 'GeneName','Count',  'Total', 'Frequency', 'Transcript_Name']]
    result2 = result2.astype({'Count': int, 'Total': int})

    temp_dict = result2.to_dict(orient='record')

    return render_template('cspec_cMDTs.html', temp_dict=temp_dict)
#
#@app.route("/fre_cMDTs")
#def fre_cMDTs():
#
#    cur = mysql.get_db().cursor()
#    cur.execute(''' SELECT cMDT, Frequency, CancerType, mart_export.`Associated Transcript Name` FROM tables4 LEFT JOIN mart_export ON tables4.cMDT = mart_export.`Ensembl Transcript ID` ORDER BY Frequency DESC''')
#    TableS4_tuple = cur.fetchall()
#    result = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','CancerType', 'Transcript_Name'])
#    temp_dict = result.to_dict(orient='record')
#
#
#    return render_template('fre_cMDTs.html', temp_dict=temp_dict)
#

@app.route("/download", methods=['GET', 'POST'])
def download():

    return render_template("download.html")

@app.route("/help", methods=['GET', 'POST'])
def help():

    cur = mysql.get_db().cursor()
    cur.execute(''' SELECT CancerType, Total FROM tables4 ''')
    TableS4_tuple = cur.fetchall()
    TableS4 = pd.DataFrame(TableS4_tuple, columns=['CancerType','Total'])
    TableS4 =  TableS4.rename(columns={'CancerType': 'CancerType','Total':'Total'})
    TableS4_uniq = TableS4.drop_duplicates()

    graphJSON =  sample_size(TableS4_uniq)

    return render_template("help.html", TableS4=TableS4, graphJSON=graphJSON)

def sample_size(TableS4):

    x = TableS4.sort_values(by=['Total']).iloc[:,0].str.split(".", n=1, expand=True).iloc[:,1]
    y = TableS4.sort_values(by=['Total']).iloc[:,1]
    data = [go.Bar(x=x, y=y, marker_color = 'indianred')]
    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route("/Cancer", methods=['GET', 'POST'])
def Cancer():

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, GeneName1, CancerSampleId, DomCancerTrans, mart_export.`Associated Transcript Name` FROM interactiondisruptionindominanttranscripts
    LEFT JOIN mart_export ON interactiondisruptionindominanttranscripts.DomCancerTrans = mart_export.`Ensembl Transcript ID`
    WHERE interactiondisruptionindominanttranscripts.Tissue = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'Associated Transcript Name'])
    df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
    df_iso = df[['CancerType2', 'Tissue', 'CancerSampleId', 'GeneName1', 'DomCancerTrans', 'Associated Transcript Name']]
    df_iso2 = df_iso.drop_duplicates()
    result =  df_iso2.rename(columns={'Associated Transcript Name':'Transcript_Name'})
    temp_dict = result.to_dict(orient='records')

    ## second table in the page representing second graph
    cur2 = mysql.get_db().cursor()
    sql2 = '''SELECT CancerType, SampleID, NumberOfMDTs FROM mdts_vs_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    dataset_muts = cur2.fetchall()

    muts_df = pd.DataFrame(dataset_muts, columns= ['CancerType', 'SampleID', 'NumberOfMDTs' ])
    muts_df[['Splitted','CancerType2']] = muts_df.CancerType.str.split('.', expand=True)

    dataset = muts_df.to_dict(orient='records')
    return render_template("Cancer_Based.html", SampleCancerType = SampleCancerType, data = temp_dict, data2 = dataset)

def CancerSpecific(SampleCancerType):

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()
    sql = ''' SELECT cMDT, Frequency, CancerType, mart_export.`Associated Transcript Name` FROM tables4 LEFT JOIN mart_export ON tables4.cMDT = mart_export.`Ensembl Transcript ID`
            WHERE tables4.CancerType = %s ORDER BY Frequency DESC LIMIT 10 '''
    adr = (SampleCancerType, )
    cur.execute(sql,adr)
    TableS4_tuple = cur.fetchall()
    result = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','CancerType', 'Transcript_Name'])
    #result = TableS4.rename(columns={'Associated Transcript Name':'Transcript_Name'})
    #result = TableS4_new.sort_values(by=['Frequency'],ascending=False).iloc[0:10,:]

    cur2 = mysql.get_db().cursor()
    sql2 = ''' SELECT CancerType, SampleID, NumberOfMDTs  FROM mdts_vs_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    cMDT_mdts_muts = cur2.fetchall()
    cMDT_dist = pd.DataFrame(cMDT_mdts_muts, columns=['CancerType', 'SampleID','NumberOfMDTs'])

    plot = make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3], subplot_titles=("Top 10 Transcripts", "Distribiton of cMDTs Across Samples"))

    trace1 = go.Bar(
                name = '% of ENSTs across ' + SampleCancerType,
                marker_color = 'indianred',
                x=result.iloc[:,3],
                y=result.iloc[:,1]*100
            )

    trace2 = go.Box(y=cMDT_dist.NumberOfMDTs,boxpoints='outliers',jitter=0.3, name = 'Distribution of cMDT Counts in Samples',
                    marker_color = '#00CC96')


    plot.append_trace(trace1, row=1, col=1)
    plot.update_yaxes(title_text="Occurence of Transcripts in Samples (%)", row=1, col=1)

    plot.append_trace(trace2, row=1, col=2)
    plot.update_yaxes(title_text="cMDT Counts across samples", row=1, col=2)
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

    genename = request.args.get('gene')
    enstid = request.args.get('enst')

    # Query interactiondisruptionindominanttranscripts from mysql table
    cur2 = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactiondisruptionindominanttranscripts
                LEFT JOIN mart_export ON interactiondisruptionindominanttranscripts.DomCancerTrans = mart_export.`Ensembl Transcript ID`
                WHERE interactiondisruptionindominanttranscripts.DomCancerTrans = %s '''
    adr = (enstid,)
    cur2.execute(sql, adr)
    isonet_tuple = cur2.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2',
                                            'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans',
                                            'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans',
                                            'StringDensityRank1','Region1', 'Transcript_Name'])

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
        # Query cancer gene census genes from mysql table
        cur4 = mysql.get_db().cursor()
        cur4.execute( '''SELECT `Gene Symbol` FROM cancer_gene_census ''' )
        cancer_gene_census_tuple = cur4.fetchall()
        df_cgc = pd.DataFrame(cancer_gene_census_tuple, columns=['Gene Symbol'])
        df_cgc =  df_cgc.rename(columns={'Gene Symbol': 'GeneName'})

        df_cgc_list = df_cgc['GeneName'].tolist()

        df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
        df = df.drop_duplicates()
        data_dict = df.to_dict(orient='records')

        #make a table for some statistics
        statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans']].drop_duplicates()
        statistics_table_dict = statistic_table.to_dict(orient='records')

        string_score = statistic_table.iloc[0,6]
        #string_score = float("{:.2f}".format(string_score))*100
        string_score = int(statistic_table.iloc[0,6]*100)
        ### DRAW PIE CHARTS

        data = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]],
                            subplot_titles=("% Interaction Lost", "Cancer Types"))


        data.add_trace(go.Pie(labels=df.drop_duplicates(subset=['CancerSampleId', 'Tissue']).Tissue), 1, 2)

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


        ## end of PIE CHART
        ## Take missed interactions for the interested ENST id from the missed_interactions table from mysqldb
        cur3 = mysql.get_db().cursor()
        sql3 = ''' SELECT ENSG, ENSP, ENST, MissInts FROM missed_interactions WHERE ENST = %s '''
        adr3 = (enstid,)
        cur3.execute(sql3, adr3)
        missed_interactions_tuple = cur3.fetchall()
        missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'MissInts'])
        Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
        Isoform_Int_Network = pd.concat([missed_interactions, Isoform_Int_Network_splitted], axis=1)

        ENSP_list = list()
        ensp_frame = list()

        for eachcolumn in range(3, len(Isoform_Int_Network.iloc[0,:])):
            Isoform_Int_Network.iloc[0,eachcolumn] = str(Isoform_Int_Network.iloc[0,eachcolumn])

            if "ENSP" in Isoform_Int_Network.iloc[0,eachcolumn]:
                ENSP_list.append(Isoform_Int_Network.iloc[0,eachcolumn])

            else:
                continue

            ensp_frame.append(ENSP_list)

        ensp_frame = functools.reduce(operator.iconcat, ensp_frame, [])

        partner_genenames = []
        #cgc_partners = []
        #cur_ensp = mysql.get_db().cursor()
        #for eachpartner in ensp_frame:
        #    try:
        #        cur_ensp_sql = '''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid = %s '''
        #        adr_ensp = (eachpartner,)
        #        cur_ensp.execute(cur_ensp_sql, adr_ensp)
        #        ensp_tuple = cur_ensp.fetchall()
        #        ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
        #        partner_genenames.append(list(ensp_genename['GeneName']))
        #        #genenames = ensp_genename_dict['GeneName'][eachpartner]
        #        #partner_genenames.append(genenames)
        #    except:
        #        pass

        try:
            placeholders = ','.join(['%s'] * len(ensp_frame))
            cur_ensp = mysql.get_db().cursor()
            cur_ensp.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des WHERE ENSPid IN (%s)'''%placeholders, tuple(ensp_frame))
            ensp_tuple = cur_ensp.fetchall()
            ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])
            partner_genenames = list(ensp_genename['GeneName'])
        except:
            pass

        cgc_partners = [eachpartner for eachpartner in partner_genenames if eachpartner in df_cgc_list]

        #for eachpartner in partner_genenames:
        #    if eachpartner in df_cgc_list:
        #        cgc_partners.append(eachpartner)

        cgc_partners_df = pd.DataFrame({'GeneName':cgc_partners})
        df_cgc_dict  = cgc_partners_df.drop_duplicates().to_dict(orient='records')

        #return render_template('network_lochen.html', genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, result=result, cgc=df_cgc_dict, df_cgc_all_dict=df_cgc_all_dict, graphJSON2=graphJSON2)
    return render_template('network.html', string_score=string_score, transcript_name=transcript_name, genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, cgc=df_cgc_dict, graphJSON2=graphJSON2)


    ### Gene Based

@app.route("/Gene", methods=['GET', 'POST'])
def Gene():

    # Take genename and enstid from url
    genename = request.args.get('gene')

    # Query cancer gene census genes from mysql table
    cur4 = mysql.get_db().cursor()
    cur4.execute( '''SELECT `Gene Symbol` FROM cancer_gene_census ''' )
    cancer_gene_census_tuple = cur4.fetchall()
    df_cgc = pd.DataFrame(cancer_gene_census_tuple, columns=['Gene Symbol'])
    df_cgc =  df_cgc.rename(columns={'Gene Symbol': 'GeneName'})
    df_cgc_list = df_cgc['GeneName'].tolist()

    cur5 = mysql.get_db().cursor()
    sql5 = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1, mart_export.`Associated Transcript Name` FROM interactiondisruptionindominanttranscripts
    LEFT JOIN mart_export ON interactiondisruptionindominanttranscripts.DomCancerTrans = mart_export.`Ensembl Transcript ID`
    WHERE interactiondisruptionindominanttranscripts.GeneName1 = %s '''
    adr5 = (genename,)
    cur5.execute(sql5, adr5)
    isonet_tuple = cur5.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1', 'Region1', 'Transcript_Name'])
    df =  df.rename(columns={'ENSG': 'ENSGid', 'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions'})
    df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
    df = df.drop_duplicates()

    df['Domain1'].replace({"-":"None"},inplace=True)
    df['Domain2'].replace({"-":"None"},inplace=True)
    df['MissedInteractions'].replace({-1:0},inplace=True)

    result = df.to_dict(orient='records')

    statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans', 'Transcript_Name']].drop_duplicates()
    statistics_table_dict = statistic_table.to_dict(orient='records')

    data_dict = df.to_dict(orient='records')

    string_score = statistic_table.iloc[0,6]
    string_score = float("{:.2f}".format(string_score))*100
    ### DRAW PIE CHART

    data = make_subplots(rows=1, cols=1, specs=[[{'type':'domain'}]])

    data.add_trace(go.Pie(labels=df.Tissue, title="Cancer Types"), 1, 1)
    data.update_traces(selector=dict(type='pie'), textinfo='label+text',hoverinfo='label+percent',
                              marker=dict(line=dict(color='#000000', width=4)))

    data.update_layout(title=" Gene Name: {} ".format(genename), showlegend=False, title_font_size=24)
    graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('GeneBased.html', string_score=string_score, df_cgc_list=df_cgc_list, genename=genename, data=data_dict, data_statistics = statistics_table_dict, result=result, graphJSON2=graphJSON2)



@app.route("/Sample", methods=['GET', 'POST'])
def Sample():

    CancerSampleId = request.args.get('sampleid')
    genename = request.args.get('gene')
    tissue = request.args.get('tissue')

    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, ENSG, GeneName1, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactiondisruptionindominanttranscripts WHERE CancerSampleId = %s '''
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
            else:
                continue

    normal_trans_id_dict = dict()
    for index,value in enumerate(normal_trans_id_list):
      normal_trans_id_dict[index] = value

    return render_template("Sample_Based.html", CancerSampleId=CancerSampleId, tissue=tissue, genename=genename, normal_trans_id_dict=normal_trans_id_dict, normal_trans_id_list=normal_trans_id_list, cancer_trans_id=cancer_trans_id)

def update_fig(CancerSampleId, genename, tissue):
    tissue = request.args.get('tissuetype')
    tissuetype = tissue.split('.')[1].replace('-','_')

    CancerSampleId = request.args.get('CanSampleId')
    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, ENSG, GeneName1, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactiondisruptionindominanttranscripts WHERE CancerSampleId = %s '''
    adr = (CancerSampleId, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()

    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
    df =  df.rename(columns={'ENSG': 'ENSGid'})

    dff = df[df.GeneName1 == genename]
    #col_name = dff.iloc[0,3] ## sample id
    cancer_trans_id = dff.iloc[0,4] ## DomCancerTrans id
    normal_trans_ids = dff['GTExMDIs'].str.split(";", expand = True) # n=1
    normal_trans_id_list = []

     #list of normal trans id

    for i in range(len(normal_trans_ids.iloc[0,:])):
        for j in range(len(normal_trans_ids.iloc[:,0])):
            if "ENST" in str(normal_trans_ids.iloc[j,i]):
                normal_trans_id_list.append(normal_trans_ids.iloc[j,i].split(":", 1)[0])
            else:
                continue

#    normal_trans_id = normal_trans_ids.iloc[0,0] # Normal Transcript Id

    #list of normal trans id

    #gtex = '_gtex'
    #my_tissue = tissuetype+gtex
    #gtex_cur_normal = mysql.get_db().cursor()
    #gtex_sql_normal = ''' SELECT * FROM ''' + my_tissue + ''' WHERE Feature REGEXP %s '''
    #gtex_adr_normal = (normal_trans_id,)
    #gtex_cur_normal.execute(gtex_sql_normal, gtex_adr_normal)
    #df_gtex_tuple = gtex_cur_normal.fetchall()
    #df_gtex_normal = pd.DataFrame(df_gtex_tuple)

    # extract normal transcripts expressions from gtex data - it can be more than 1 transcript
    gtex = '_gtex'
    my_tissue = tissuetype+gtex
    #placeholders = '|'.join(['%s'] * len(normal_trans_id_list))
    gtex_cur_normal = mysql.get_db().cursor()
    df_gtex_normal = pd.DataFrame()
    for eachtranscript in normal_trans_id_list:

        gtex_sql_normal = '''SELECT * FROM ''' + my_tissue + ''' WHERE Feature REGEXP %s'''
        gtex_adr_normal = (eachtranscript,)
        gtex_cur_normal.execute(gtex_sql_normal, gtex_adr_normal)
        df_gtex_tuple = gtex_cur_normal.fetchall()
        df_gtex_normal_df = pd.DataFrame(df_gtex_tuple)
        df_gtex_normal = df_gtex_normal.append(df_gtex_normal_df, ignore_index = True)
    #df_gtex_normal_exp_list = list(df_gtex_normal['Feature'])

    # extract cancer transcript expression from gtex data - only 1 transcript

    gtex_cur_cancer = mysql.get_db().cursor()
    gtex_sql_cancer = ''' SELECT * FROM ''' + my_tissue + ''' WHERE Feature REGEXP %s '''
    gtex_adr_cancer = (cancer_trans_id,)
    gtex_cur_cancer.execute(gtex_sql_cancer, gtex_adr_cancer)
    df_gtex_tuple_2 = gtex_cur_cancer.fetchall()
    df_gtex_cancer = pd.DataFrame(df_gtex_tuple_2)

    pcawg = '_pcawg'
    my_tissue2 = tissuetype + pcawg

    #pcawg_cur_normal = mysql.get_db().cursor()
    #pcawg_sql_normal = ''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM ''' + my_tissue2 + ''' WHERE Feature REGEXP %s '''
    #pcawg_adr_normal = (normal_trans_id,)
    #pcawg_cur_normal.execute(pcawg_sql_normal, pcawg_adr_normal)
    #df_pcawg_tuple = pcawg_cur_normal.fetchall()
    #df_pcawg_normal = pd.DataFrame(df_pcawg_tuple)

    #extract normal transcripts expressions from pcawg data - it can be more than 1 transcript
    #placeholders = '|'.join(['%s'] * len(normal_trans_id_list))
    pcawg_cur_normal = mysql.get_db().cursor()
    df_pcawg_normal = pd.DataFrame()
    for eachtranscript in normal_trans_id_list:
        pcawg_sql_normal = ''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM ''' + my_tissue2 + ''' WHERE Feature REGEXP %s '''
        pcawg_adr_normal = (eachtranscript,)
        pcawg_cur_normal.execute(pcawg_sql_normal, pcawg_adr_normal)
        df_pcawg_tuple = pcawg_cur_normal.fetchall()
        df_pcawg_normal_df = pd.DataFrame(df_pcawg_tuple)
        df_pcawg_normal = df_pcawg_normal.append(df_pcawg_normal_df, ignore_index = True)
       # pcawg_cur_normal.execute('''SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM ''' + my_tissue2 + ''' WHERE Feature REGEXP %s'''%placeholders, normal_trans_id_list)


    # extract cancer transcript expression from pcawg data - only 1 transcript
    pcawg_cur_cancer = mysql.get_db().cursor()
    pcawg_sql_cancer = ''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM ''' + my_tissue2 + ''' WHERE Feature REGEXP %s '''
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
    #count_df2 = pd.DataFrame()

    #for each_normal in list(df_gtex_normal.iloc[:,0]):
    for each_normal in set(list(df_gtex_normal.iloc[:,0])):
        x = [cancer_trans_id]
        x.append(each_normal)

        count_data_normal_trans_median = np.median(df_gtex_normal[df_gtex_normal.iloc[:,0] == each_normal].iloc[0,1:len(df_gtex_normal.columns)])
        std_gtex_normal = statistics.pstdev(df_gtex_normal[df_gtex_normal.iloc[:,0] == each_normal].iloc[0,1:len(df_gtex_normal.columns)])

        count_df = count_df.append(pd.DataFrame({"ENST": each_normal,"GTEx":count_data_normal_trans_median, "PCAWG":float(df_pcawg_normal[df_pcawg_normal.iloc[:,0] == each_normal].iloc[0,1]), "Std":std_gtex_normal}, index=[each_normal]))

    cancer_exp = pd.DataFrame({"ENST": cancer_trans_id, "GTEx":count_data_cancer_trans_median, "PCAWG":float(df_pcawg_cancer.iloc[:,1]), "Std":std_gtex_cancer}, index=[cancer_trans_id])
    count_df = count_df.append(cancer_exp)

    for col in count_df.iloc[:,1:2].columns:
        data.add_trace(go.Bar(x=count_df.index, y=count_df[col], name = col, error_y=dict(type='data', array=list(count_df['Std'])),marker_color='indianred'))

    for col in count_df.iloc[:,2:3].columns:
        data.add_trace(go.Bar(x=count_df.index, y=count_df[col], name = col, marker_color='lightsalmon'))

    data.update_layout(title=genename, yaxis_title="TPM Count", title_font_size=24)

    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/Gene_Based2', methods=['GET', 'POST'])
def change_features4():

    CancerSampleId = request.args.get('CanSampleId')
    genename = request.args.get('genename')
    tissue = request.args.get('tissuetype')

    graphJSON = update_fig(CancerSampleId, genename, tissue)

    return graphJSON


@app.errorhandler(500)
def page_not_found(e):
    return render_template('500.html')

if __name__ == "__main__":
  app.run(debug=True)