from flask import Flask, render_template, request
from flask_bootstrap import Bootstrap
import plotly
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import json
from flask_table import Table, Col
from flaskext.mysql import MySQL
import MySQLdb.cursors
import yaml
import csv
import pymysql
import ast
from ipywidgets import widgets ## to read dictionary from a file

app = Flask(__name__)


db = yaml.load(open('db.yaml'), Loader=yaml.FullLoader)
app.config['MYSQL_DATABASE_HOST'] = db['mysql_host']
app.config['MYSQL_DATABASE_USER'] = db['mysql_user']
app.config['MYSQL_DATABASE_PASSWORD'] = db['mysql_password']
app.config['MYSQL_DATABASE_DB'] = db['mysql_db']
app.config['MYSQL_DATABASE_PORT'] = db['mysql_port']


mysql = MySQL()
mysql.init_app(app)


ensp_genename = pd.read_table('/Users/tulaykarakulak/Documents/PhD/Projects/CanIsoNet_Web/Flask/static/data/ensp_genename.txt')


@app.route("/")
def home():
    cur = mysql.get_db().cursor()
    cur.execute('''SELECT * FROM mdts_vs_muts''')
    cMDT_mdts_muts = cur.fetchall()
    cMDT_dist = pd.DataFrame(cMDT_mdts_muts, columns=['CancerType', 'SampleID','NumberOfMDTs','NumberOfMutations'])

    cancertypes = cMDT_dist.CancerType.unique()
    cancertypes2 = pd.DataFrame(cancertypes,columns=['CancerType'])
    cancertypes2[['Splitted','CancerType2']] = cancertypes2.CancerType.str.split('.', expand=True)


    cur5 = mysql.get_db().cursor()
    cur5.execute('''SELECT * FROM cancertypes''')
    cancertypes_sql = cur5.fetchall()
    cancer_types = pd.DataFrame(cancertypes_sql, columns=['PCAWG-Code', 'Cancer Type Name'])
    cancer_types = cancer_types.rename(columns={'PCAWG-Code':'CancerType2', 'Cancer Type Name':'Cancer Type Name'})

    result = pd.merge(cancer_types, cancertypes2, how='inner', on='CancerType2')

    cancertypes_dict = result.to_dict(orient='records')



    cur2 = mysql.get_db().cursor()
    cur2.execute(''' SELECT GeneName1, DomCancerTrans FROM interactiondisruptionindominanttranscripts ''')
    isonet_tuple = cur2.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['GeneName1','DomCancerTrans'])
    df = df.drop_duplicates()
    df =  df.rename(columns={'GeneName1': 'GeneName1','DomCancerTrans':'DomCancerTrans'})


    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur3.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})

    result = pd.merge(df, df_biomart, how='inner', on='DomCancerTrans')

    temp_dict = result.to_dict(orient='records')

    return render_template('home.html', data=cancertypes_dict, data2=temp_dict)


@app.route("/Cancer_Based", methods=['GET', 'POST'])
def Cancer_Based():

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()

    sql = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans FROM interactiondisruptionindominanttranscripts WHERE Tissue = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans'])
    df =  df.rename(columns={'Tissue': 'Tissue', 'ENSG': 'ENSGid','NumberOfGtexMDIs': 'NumberOfGtexMDIs','GeneName1': 'GeneName1','GeneName2': 'GeneName2',
                                'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions',
                                'Pfam1': 'Pfam1','Domain1': 'Domain1','Pfam2': 'Pfam2','Domain2': 'Domain2', 'CancerSampleId': 'CancerSampleId', 'DomCancerTrans':'DomCancerTrans'})

    df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)

    df_iso = df[['CancerType2', 'Tissue', 'CancerSampleId', 'GeneName1', 'DomCancerTrans']]
    df_iso2 = df_iso.drop_duplicates()

    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur3.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})
    result = pd.merge(df_iso2, df_biomart, how='inner', on='DomCancerTrans')
    temp_dict = result.to_dict(orient='records')

    ## second table in the page representing second graph
    cur2 = mysql.get_db().cursor()
    sql2 = '''SELECT CancerType, SampleID, NumberOfMDTs, NumberofMutations FROM mdts_vs_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    dataset_muts = cur2.fetchall()

    muts_df = pd.DataFrame(dataset_muts, columns= ['CancerType', 'SampleID', 'NumberOfMDTs' , 'NumberOfMutations'])
    muts_df[['Splitted','CancerType2']] = muts_df.CancerType.str.split('.', expand=True)

    dataset = muts_df.to_dict(orient='records')


    return render_template("Cancer_Based.html", SampleCancerType = SampleCancerType, data = temp_dict, data2 = dataset)

def CancerSpecific(SampleCancerType):

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()
    sql = ''' SELECT * FROM TableS4 WHERE CancerType = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    TableS4_tuple = cur.fetchall()
    TableS4 = pd.DataFrame(TableS4_tuple, columns=['CancerType', 'ENSG','cMDT','GeneName','Count','Total','Frequency'])

    TableS4 = TableS4.sort_values(by=['Frequency'],ascending=False).iloc[0:10,:]

    cur2 = mysql.get_db().cursor()
    sql2 = ''' SELECT * FROM mdts_vs_muts WHERE CancerType = %s '''
    adr2 = (SampleCancerType, )
    cur2.execute(sql2, adr2)
    cMDT_mdts_muts = cur2.fetchall()
    cMDT_dist = pd.DataFrame(cMDT_mdts_muts, columns=['CancerType', 'SampleID','NumberOfMDTs','NumberOfMutations'])

    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur3.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'cMDT','Associated Transcript Name':'Transcript_Name'})
    result = pd.merge(TableS4, df_biomart, how='inner', on='cMDT')

    plot = make_subplots(rows=1, cols=2, column_widths=[0.7, 0.3], subplot_titles=("Top 10 Transcripts", "Distribiton of cMDTs Across Samples"))

    trace1 = go.Bar(
                name = '% of ENSTs across ' + SampleCancerType,
                marker_color = 'indianred',
                x=result.iloc[:,7],
                y=result.iloc[:,6]*100
            )

    trace2 = go.Box(y=cMDT_dist.NumberOfMDTs,boxpoints='all',jitter=0.3, name = 'Distribution of cMDT Counts in Samples',
                    marker_color = '#00CC96')


    plot.append_trace(trace1, row=1, col=1)
    plot.append_trace(trace2, row=1, col=2)

    graphJSON = json.dumps(plot, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/CancerSpecific', methods=['GET', 'POST'])
def change_features5():

    SampleCancerType = request.args.get('cancer')
    graphJSON = CancerSpecific(SampleCancerType)

    return graphJSON

@app.route("/Gene_Based", methods=['GET', 'POST'])
def Gene_Based():

    genename = request.args.get('gene')
    enstid = request.args.get('enst')
    partner_genenames = []
    cgc_partners = []

    cur2 = mysql.get_db().cursor()
    sql = '''SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1 FROM interactiondisruptionindominanttranscripts WHERE DomCancerTrans = %s '''
    adr = (enstid,)
    cur2.execute(sql, adr)
    isonet_tuple = cur2.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1','Region1'])
    df =  df.rename(columns={'Tissue': 'Tissue', 'ENSG': 'ENSGid','NumberOfGtexMDIs': 'NumberOfGtexMDIs','GeneName1': 'GeneName1','GeneName2': 'GeneName2',
                                'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions',
                                'Pfam1': 'Pfam1','Domain1': 'Domain1','Pfam2': 'Pfam2','Domain2': 'Domain2', 'CancerSampleId': 'CancerSampleId', 'DomCancerTrans':'DomCancerTrans', 'StringDensityRank1':'StringDensityRank1', 'Region1':'Region1'})


    cur4 = mysql.get_db().cursor()
    cur4.execute( '''SELECT `Gene Symbol` FROM cancer_gene_census ''' )
    cancer_gene_census_tuple = cur4.fetchall()
    df_cgc = pd.DataFrame(cancer_gene_census_tuple, columns=['Gene Symbol'])
    df_cgc =  df_cgc.rename(columns={'Gene Symbol': 'GeneName'})
    df_cgc_all_dict = df_cgc.to_dict(orient='records')

    df_cgc_list = df_cgc['GeneName'].tolist()


    ## Take missed interactions for the ENST id from the missed_interactions table from mysql database
    if enstid != None:

        df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
        df = df.drop_duplicates()

        statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1']].drop_duplicates()
        statistics_table_dict = statistic_table.to_dict(orient='records')


        ### DRAW PIE CHART
        data = make_subplots(rows=3, cols=1, specs=[[{"type": "pie"}],[{"type": "pie"}],[{"type": "pie"}]])

        data.add_trace(go.Pie(values=[statistic_table.iloc[0,6]*100,(100-(statistic_table.iloc[0,6]*100))], title="String Density Score"), 1, 1)
        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(colors=['lightgreen', 'darkorange'], line=dict(color='#000000', width=4)))
        data.add_trace(go.Pie(values=[(statistic_table.iloc[0,2]-statistic_table.iloc[0,3])*100/48,statistic_table.iloc[0,3]*100/48], labels=['# of Remaining Interaction', '# of  Interaction Lost'], title="Interaction Lost"),2,1)
        data.update_traces(selector=dict(type='pie'),
                          marker=dict(colors=['gold', 'mediumturquoise'], line=dict(color='#000000', width=4)))

        data.add_trace(go.Pie(labels=df.Tissue, title="Cancer Types"), 3, 1)
        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(line=dict(color='#000000', width=4)))

        data.update_layout(showlegend=False, title_font_size=5)
        graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

        ## end of PIE CHART

        data_dict = df.to_dict(orient='records')


## to put the transcript names inside the table
        cur6 = mysql.get_db().cursor()
        cur6.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
        biomart_tuple = cur6.fetchall()
        df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
        df_biomart = df_biomart.drop_duplicates()
        df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})

        result_df = pd.merge(df, df_biomart, how='inner', on='DomCancerTrans')
        result = result_df.to_dict(orient='records')

## end of the table


        cur3 = mysql.get_db().cursor()
        sql3 = ''' SELECT ENSG, ENSP, ENST, MissInts FROM missed_interactions WHERE ENST = %s '''
        adr3 = (enstid,)
        cur3.execute(sql3, adr3)
        missed_interactions_tuple = cur3.fetchall()
        missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'MissInts'])
        missed_interactions =  missed_interactions.rename(columns={'ENSG':'ENSG', 'ENSP':'ENSP', 'ENST':'ENST', 'MissInts': 'MissInts'})

        Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
        frames = [missed_interactions, Isoform_Int_Network_splitted]
        Isoform_Int_Network = pd.concat(frames, axis=1)
        ## Make a dictionary where ENSG ids are key that has ENST ids as first value in each list and missed interactions for each ENST id are the value

        ENSP_list = list()
        ensp_frame = list()
        for eachcolumn in range(3, len(Isoform_Int_Network.iloc[0,:])):
            Isoform_Int_Network.iloc[0,eachcolumn] = str(Isoform_Int_Network.iloc[0,eachcolumn])

            if "ENSP" in Isoform_Int_Network.iloc[0,eachcolumn]:
                ENSP_list.append(Isoform_Int_Network.iloc[0,eachcolumn])

            else:
                continue

            ensp_frame.append(ENSP_list)

        subset_missed_int2 = pd.DataFrame({'Miss_Ints': ensp_frame })
        subset_missed_int3 = pd.concat([Isoform_Int_Network,subset_missed_int2], axis=1)

        dictionary_missed_ints = subset_missed_int3.groupby('ENSG')[['ENST','Miss_Ints']].apply(lambda g: list(map(tuple, g.values.tolist()))).to_dict()

        for eachkey in dictionary_missed_ints:

            for i in range(0,len(dictionary_missed_ints[eachkey])):

                if enstid == dictionary_missed_ints[eachkey][i][0]:
                    partner_enspids = [dictionary_missed_ints[eachkey][i][1]]
                    make_df = pd.DataFrame(partner_enspids)
                    for i in range(0, len(make_df.iloc[0,:])):
                            for j in range(0, len(ensp_genename)):
                                if make_df.iloc[0,i] == ensp_genename.iloc[j,0]:
                                        partner_genenames.append(ensp_genename.iloc[j,1])

    else:
        cur5 = mysql.get_db().cursor()
        sql5 = '''SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1 FROM interactiondisruptionindominanttranscripts WHERE GeneName1 = %s '''
        adr5 = (genename,)
        cur5.execute(sql5, adr5)
        isonet_tuple = cur5.fetchall()
        df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'StringDensityRank1', 'Region1'])
        df =  df.rename(columns={'Tissue': 'Tissue', 'ENSG': 'ENSGid','NumberOfGtexMDIs': 'NumberOfGtexMDIs','GeneName1': 'GeneName1','GeneName2': 'GeneName2',
                                    'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions',
                                    'Pfam1': 'Pfam1','Domain1': 'Domain1','Pfam2': 'Pfam2','Domain2': 'Domain2', 'CancerSampleId': 'CancerSampleId', 'DomCancerTrans':'DomCancerTrans', 'StringDensityRank1':'StringDensityRank1', 'Region1':'Region1'})


        df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
        df = df.drop_duplicates()

        data_dict = df.to_dict(orient='records')

## to put the transcript names inside the table
        cur6 = mysql.get_db().cursor()
        cur6.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
        biomart_tuple = cur6.fetchall()
        df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
        df_biomart = df_biomart.drop_duplicates()
        df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})

        result_df = pd.merge(df, df_biomart, how='inner', on='DomCancerTrans')
        result = result_df.to_dict(orient='records')
## end of the table


        statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1']].drop_duplicates()
        statistics_table_dict = statistic_table.to_dict(orient='records')

        data_dict = df.to_dict(orient='records')


        ### DRAW PIE CHART
        #data = px.pie(values=[statistic_table.iloc[0,6]*100,(100-(statistic_table.iloc[0,6]*100))])
        data = make_subplots(rows=3, cols=1, specs=[[{"type": "pie"}],[{"type": "pie"}],[{"type": "pie"}]])

        data.add_trace(go.Pie(values=[statistic_table.iloc[0,6]*100,(100-(statistic_table.iloc[0,6]*100))], title="String Density Score"), 1, 1)
        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(colors=['lightgreen', 'darkorange'],
                                  line=dict(color='#000000', width=4)))
        data.add_trace(go.Pie(values=[(statistic_table.iloc[0,2]-statistic_table.iloc[0,3])*100/48,
                        statistic_table.iloc[0,3]*100/48],
                        labels=['# of Remaining Interaction', '# of  Interaction Lost'],
                        title="Interaction Lost"),2,1)
        data.update_traces(selector=dict(type='pie'),
                          marker=dict(colors=['gold', 'mediumturquoise'],
                          line=dict(color='#000000', width=4)))

        data.add_trace(go.Pie(labels=df.Tissue, title="Cancer Types"), 3, 1)
        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(line=dict(color='#000000', width=4)))

        data.update_layout(showlegend=False, title_font_size=5)
        graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    for eachpartner in partner_genenames:
        if eachpartner in df_cgc_list:
            cgc_partners.append(eachpartner)

    cgc_partners_df = pd.DataFrame({'GeneName':cgc_partners})
    df_cgc_dict  = cgc_partners_df.to_dict(orient='records')

    return render_template('network_lochen.html', genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, result=result, cgc=df_cgc_dict, df_cgc_all_dict=df_cgc_all_dict, graphJSON2=graphJSON2)


@app.route("/help", methods=['GET', 'POST'])
def help():

    cur = mysql.get_db().cursor()
    cur.execute(''' SELECT CancerType, Total FROM TableS4 ''')
    TableS4_tuple = cur.fetchall()
    TableS4 = pd.DataFrame(TableS4_tuple, columns=['CancerType','Total'])
    TableS4 =  TableS4.rename(columns={'CancerType': 'CancerType','Total':'Total'})
    TableS4_uniq = TableS4.drop_duplicates()

    graphJSON =  sample_size(TableS4)

    return render_template("help.html", TableS4=TableS4, graphJSON=graphJSON)

def sample_size(TableS4):

    data = [go.Bar(x=TableS4.iloc[:,0], y=TableS4.iloc[:,1], marker_color = 'indianred')]
    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route("/Sample_Based", methods=['GET', 'POST'])
def Sample_Based():

    CancerSampleId = request.args.get('CanSampleId')
    genename = request.args.get('genename')
    tissue = request.args.get('tissuetype')

    return render_template("Sample_Based.html", CancerSampleId=CancerSampleId, tissue=tissue, genename=genename)

def update_fig(CancerSampleId, genename, tissue):
    tissue = request.args.get('tissuetype')
    tissuetype = tissue.split('.')[1].replace('-','_')

    CancerSampleId = request.args.get('CanSampleId')
    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactiondisruptionindominanttranscripts WHERE CancerSampleId = %s '''
    adr = (CancerSampleId, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()

    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2', 'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans', 'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
    df =  df.rename(columns={'Tissue': 'Tissue', 'ENSG': 'ENSGid','NumberOfGtexMDIs': 'NumberOfGtexMDIs','GeneName1': 'GeneName1','GeneName2': 'GeneName2',
                                'TotalNumberOfStringInt': 'NumberOfStringInt','NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions',
                                'Pfam1': 'Pfam1','Domain1': 'Domain1','Pfam2': 'Pfam2','Domain2': 'Domain2', 'CancerSampleId': 'CancerSampleId', 'DomCancerTrans':'DomCancerTrans', 'GTExMDIs':'GTExMDIs'})


    cur_gtex = mysql.get_db().cursor()
    gtex = '_gtex'
    my_tissue = tissuetype+gtex
    cur_gtex.execute( ''' SELECT * FROM ''' + my_tissue )

    gtex_tuple = cur_gtex.fetchall()

    x = '/Users/tulaykarakulak/Documents/PhD/Projects/CanIsoNet_Web/Flask_IsoNet/datasets/' + tissue.split('.')[0] + '/' + tissue.split('.')[1] + '/' + 'gtex_column_names'
    column_names = pd.read_table(x, header=None)[0].to_list()

    df_gtex = pd.DataFrame(gtex_tuple, columns = column_names)



    cur_pcawg = mysql.get_db().cursor()
    pcawg = '_pcawg'
    my_tissue2 = tissuetype + pcawg
    cur_pcawg.execute ( ''' SELECT * FROM ''' + my_tissue2)

    pcawg_tuple = cur_pcawg.fetchall()


    x = '/Users/tulaykarakulak/Documents/PhD/Projects/CanIsoNet_Web/Flask_IsoNet/datasets/' + tissue.split('.')[0] + '/' + tissue.split('.')[1] + '/' + 'pcawg_column_names'
    column_names_pcawg = pd.read_table(x, header=None)[0].to_list()

    df_pcawg = pd.DataFrame(pcawg_tuple, columns =  column_names_pcawg)


    data = go.Figure()

    dff = df[df.GeneName1 == genename]
    new_data2 = pd.DataFrame()
    new_data3 = pd.DataFrame()
    new_data = pd.DataFrame()
    for i in range(0, len(dff.iloc[:,0])):
        tissue_type = dff.iloc[0,0]
        col_name = dff.iloc[0,11] ## sample id
        cancer_trans_id = dff.iloc[i,12] ## DomCancerTrans id
        normal_trans_ids = dff['GTExMDIs'].str.split(":", n = 1, expand = True)
        normal_trans_id = normal_trans_ids.iloc[i,0] # Normal Transcript Id

                ## Transcript Counts of Transcripts in pcawg data
        count_data_cancer_trans = df_pcawg[df_pcawg.iloc[:,0].str.split(".", n = 1, expand = True)[0] == cancer_trans_id][col_name]
        count_data_normal_trans = df_pcawg[df_pcawg.iloc[:,0].str.split(".", n = 1, expand = True)[0] == normal_trans_id][col_name]

        x=[cancer_trans_id, normal_trans_id]
        y=[float(count_data_cancer_trans),float(count_data_normal_trans)]

        new_data_gtex = df_gtex[df_gtex.iloc[:,0].str.split(".", n = 1, expand = True)[0] == normal_trans_id]
        new_data_cancer = df_gtex[df_gtex.iloc[:,0].str.split(".", n = 1, expand = True)[0] == cancer_trans_id]

        count_data_cancer_trans_median = np.median(new_data_cancer.iloc[0,1:len(new_data_cancer.columns)])
        count_data_normal_trans_median = np.median(new_data_gtex.iloc[0,1:len(new_data_gtex.columns)])
        std_cancer = np.std(new_data_cancer.iloc[0,1:len(new_data_cancer.columns)])
        std_normal = np.std(new_data_gtex.iloc[0,1:len(new_data_gtex.columns)])

        x2=[cancer_trans_id, normal_trans_id]
        y2=[count_data_cancer_trans_median,count_data_normal_trans_median]

    data.add_trace(go.Bar(
        name='PCAWG',
        marker_color='indianred',
        x=x,
        y=y)
        )

    data.add_trace(go.Bar(
        name='GTEx',
        marker_color='lightsalmon',
        x=x2,
        y=y2)
        )


    graphJSON = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON


@app.route('/Gene_Based2', methods=['GET', 'POST'])
def change_features4():

    CancerSampleId = request.args.get('CanSampleId')
    genename = request.args.get('genename')
    tissue = request.args.get('tissuetype')

    graphJSON = update_fig(CancerSampleId, genename, tissue)

    return graphJSON


if __name__ == "__main__":
  app.run(debug=True)
