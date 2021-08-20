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

app = Flask(__name__)
#app = Quart(__name__)


app.config['MYSQL_DATABASE_HOST'] = 'abxka.mysql.eu.pythonanywhere-services.com'
app.config['MYSQL_DATABASE_USER'] = 'abxka'
app.config['MYSQL_DATABASE_PASSWORD'] = 'CanIsoNetMySQL?='
app.config['MYSQL_DATABASE_DB'] = 'abxka$Isonet'

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


#    cur3 = mysql.get_db().cursor()
#    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
#    biomart_tuple = cur3.fetchall()
#    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
#    df_biomart = df_biomart.drop_duplicates()
#    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})
#    result = genenamecmdt.join(df_biomart.set_index('DomCancerTrans'), on='DomCancerTrans')
#    temp_dict2 = result.to_dict(orient='records')

    return render_template('home.html', data=cancertypes_dict, data2=genenamecmdt_gene_list, data3=genenamecmdt_gene_cmdt, temp_dict=temp_dict, genenamecmdt_gene_tn=genenamecmdt_gene_tn)

@app.route("/Cancer_Based", methods=['GET', 'POST'])
def Cancer_Based():

    SampleCancerType = request.args.get('cancer')

    cur = mysql.get_db().cursor()

    #sql = ''' SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans FROM interactiondisruptionindominanttranscripts WHERE Tissue = %s '''
    sql = ''' SELECT Tissue, GeneName1, CancerSampleId, DomCancerTrans FROM interactiondisruptionindominanttranscripts WHERE Tissue = %s '''
    #sql = ''' CREATE UNIQUE INDEX TissueType ON interactiondisruptionindominanttranscripts(Tissue(500)) '''
    #cur.execute(sql)
    #sql = ''' SELECT Tissue, GeneName1, CancerSampleId, DomCancerTrans FROM interactiondisruptionindominanttranscripts USE INDEX (TissueType) WHERE Tissue = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'GeneName1', 'CancerSampleId', 'DomCancerTrans'])
    #df =  df.rename(columns={'Tissue': 'Tissue','GeneName1': 'GeneName1', 'CancerSampleId': 'CancerSampleId', 'DomCancerTrans':'DomCancerTrans'})
    df[['Splitted','CancerType2']] = df.Tissue.str.split('.', expand=True)
    df_iso = df[['CancerType2', 'Tissue', 'CancerSampleId', 'GeneName1', 'DomCancerTrans']]
    df_iso2 = df_iso.drop_duplicates()

    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur3.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})
    #result = pd.merge(df_iso2, df_biomart, how='inner', on='DomCancerTrans')
    result = df_iso2.join(df_biomart.set_index('DomCancerTrans'), on='DomCancerTrans')
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
    sql = ''' SELECT cMDT, Frequency, CancerType FROM tables4 WHERE CancerType = %s '''
    adr = (SampleCancerType, )
    cur.execute(sql, adr)
    TableS4_tuple = cur.fetchall()
    TableS4 = pd.DataFrame(TableS4_tuple, columns=['cMDT','Frequency','CancerType'])
    TableS4 = TableS4.sort_values(by=['Frequency'],ascending=False).iloc[0:10,:]

    cur3 = mysql.get_db().cursor()
    cur3.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur3.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    #df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'cMDT','Associated Transcript Name':'Transcript_Name'})
    result = pd.merge(TableS4, df_biomart, how='inner', on='cMDT')

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

    trace2 = go.Box(y=cMDT_dist.NumberOfMDTs,boxpoints='all',jitter=0.3, name = 'Distribution of cMDT Counts in Samples',
                    marker_color = '#00CC96')


    plot.append_trace(trace1, row=1, col=1)
    plot.update_yaxes(title_text="Occurence of Transcripts in Samples (%)", row=1, col=1)

    plot.append_trace(trace2, row=1, col=2)
    plot.update_yaxes(title_text="cMDT Counts across samples", row=1, col=2)

    graphJSON = json.dumps(plot, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON

@app.route('/CancerSpecific', methods=['GET', 'POST'])
def change_features5():

    SampleCancerType = request.args.get('cancer')
    graphJSON = CancerSpecific(SampleCancerType)

    return graphJSON

@app.route("/Transcript_Based", methods=['GET', 'POST'])
def Transcript_Based():

    genename = request.args.get('gene')
    enstid = request.args.get('enst')

    cur_ensp = mysql.get_db().cursor()
    #sql_ensp = '''SELECT ENSGid, ENSTid, GeneName, TranscriptName FROM ensg_enst_ensp_des WHERE ENSTid = %s '''
    cur_ensp.execute('''SELECT ENSPid, GeneName  FROM ensg_enst_ensp_des ''')
    #adr_ensp = (enstid,)
    #cur_ensp.execute(sql_ensp, adr_ensp)
    ensp_tuple = cur_ensp.fetchall()
    ensp_genename = pd.DataFrame(ensp_tuple, columns=['ENSPid', 'GeneName'])

    #ensg_enst_ensp_description = pd.read_table('/home/abxka/ensg_enst_ensp_des.txt')
    #ensp_genename = ensg_enst_ensp_description.iloc[:,2:4]

    # Take genename and enstid from url
    ensp_genenames = ensp_genename.set_index('ENSPid', drop=False)
    ensp_genename_dict = ensp_genenames.to_dict('series')
    partner_genenames = []
    cgc_partners = []

    # Query interactiondisruptionindominanttranscripts from mysql table
    cur2 = mysql.get_db().cursor()
    sql = '''SELECT Tissue, ENSG, NumberOfGtexMDIs, GeneName1, GeneName2, TotalNumberOfStringInt, NumberOfUniqMissedInteractionsOfDomCancerTrans, Pfam1, Domain1, Pfam2, Domain2, CancerSampleId, DomCancerTrans, StringDensityRank1, Region1 FROM interactiondisruptionindominanttranscripts WHERE DomCancerTrans = %s '''
    adr = (enstid,)
    cur2.execute(sql, adr)
    isonet_tuple = cur2.fetchall()
    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'NumberOfGtexMDIs', 'GeneName1', 'GeneName2',
                                            'TotalNumberOfStringInt', 'NumberOfUniqMissedInteractionsOfDomCancerTrans',
                                            'Pfam1', 'Domain1', 'Pfam2', 'Domain2', 'CancerSampleId', 'DomCancerTrans',
                                            'StringDensityRank1','Region1'])

    df = df.rename(columns={'ENSG': 'ENSGid', 'NumberOfUniqMissedInteractionsOfDomCancerTrans': 'MissedInteractions','TotalNumberOfStringInt': 'NumberOfStringInt'})

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

        ### DRAW PIE CHARTS
        data = make_subplots(rows=1, cols=3, specs=[[{'type':'domain'}, {'type':'domain'}, {'type':'domain'}]],
                            subplot_titles=("STRING Density Score", "% Interaction Lost", "Cancer Types"))

        data.add_trace(go.Pie(values=[statistic_table.iloc[0,6]*100,(100-(statistic_table.iloc[0,6]*100))],
                            labels=['STRING Density Score',''], pull=[0.2, 0]), 1, 1)
        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(colors=['lightgreen', 'White'], line=dict(color='#000000', width=4)))

        data.add_trace(go.Pie(labels=df.drop_duplicates(subset=['CancerSampleId', 'Tissue']).Tissue), 1, 3)

        data.update_traces(selector=dict(type='pie'),
                                  marker=dict(line=dict(color='#000000', width=4)))


        if statistic_table.iloc[0,2] == 0:

            data.add_trace(go.Pie(values=[100, 0], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,2)
            data.update_traces(selector=dict(type='pie'),
                              marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))

        else:
            data.add_trace(go.Pie(values=[(statistic_table.iloc[0,2]-statistic_table.iloc[0,3])*100/(statistic_table.iloc[0,2]),statistic_table.iloc[0,3]*100/(statistic_table.iloc[0,2])], labels=['% of Remaining Interaction', '% of  Interaction Lost']),1,2)
            data.update_traces(selector=dict(type='pie'),
                              marker=dict(colors=['mediumturquoise', 'gold'], line=dict(color='#000000', width=4)))

        data.update_layout(showlegend=False, title_font_size=18)
        graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)


        ## end of PIE CHART
        ## Query Ensembl Mart and use it to merge transcript IDs in interactiondistruption table
        cur6 = mysql.get_db().cursor()
        sql6 = ''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export WHERE `Ensembl Transcript ID` = %s '''
        adr6 = (enstid,)
        cur6.execute(sql6, adr6)
        biomart_tuple = cur6.fetchall()
        df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
        df_biomart = df_biomart.drop_duplicates()
        df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})
        result_df = pd.merge(df, df_biomart, how='inner', on='DomCancerTrans')
        result = result_df.to_dict(orient='records')


        ## Take missed interactions for the interested ENST id from the missed_interactions table from mysqldb
        cur3 = mysql.get_db().cursor()
        sql3 = ''' SELECT ENSG, ENSP, ENST, MissInts FROM missed_interactions WHERE ENST = %s '''
        adr3 = (enstid,)
        cur3.execute(sql3, adr3)
        missed_interactions_tuple = cur3.fetchall()
        missed_interactions = pd.DataFrame(missed_interactions_tuple, columns=['ENSG', 'ENSP', 'ENST', 'MissInts'])
        #missed_interactions =  missed_interactions.rename(columns={'ENSG':'ENSG', 'ENSP':'ENSP', 'ENST':'ENST', 'MissInts': 'MissInts'})

        Isoform_Int_Network_splitted = pd.DataFrame(missed_interactions.MissInts.str.split(':').tolist())
        #frames = [missed_interactions, Isoform_Int_Network_splitted]
        #Isoform_Int_Network = pd.concat(frames, axis=1)
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

        #for eachpartner in ensp_frame:
        #    for j in range(0, len(ensp_genename.iloc[:,0])):
        #        if eachpartner == ensp_genename.iloc[j,0]:
        #            partner_genenames.append(ensp_genename.iloc[j,1])

        #ensp_genenames2 = ensp_genename.dropna()
        #for eachpartner in ensp_frame:
        #    try:
        #        genenames = ensp_genenames2[ensp_genenames2.iloc[:,0] == eachpartner].iloc[0,1]
        #        partner_genenames.append(genenames)
        #    except:
        #        pass


        for eachpartner in ensp_frame:
            try:
                genenames = ensp_genename_dict['GeneName'][eachpartner]
                partner_genenames.append(genenames)
            except:
                pass


        for eachpartner in partner_genenames:
            if eachpartner in df_cgc_list:
                cgc_partners.append(eachpartner)

        cgc_partners_df = pd.DataFrame({'GeneName':cgc_partners})
        df_cgc_dict  = cgc_partners_df.drop_duplicates().to_dict(orient='records')

        #return render_template('network_lochen.html', genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, result=result, cgc=df_cgc_dict, df_cgc_all_dict=df_cgc_all_dict, graphJSON2=graphJSON2)
    return render_template('network_lochen.html', genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, result=result, cgc=df_cgc_dict, graphJSON2=graphJSON2)


    ### Gene Based

@app.route("/GeneBased", methods=['GET', 'POST'])
def GeneBased():

    # Take genename and enstid from url
    genename = request.args.get('gene')

    # Query cancer gene census genes from mysql table
    cur4 = mysql.get_db().cursor()
    cur4.execute( '''SELECT `Gene Symbol` FROM cancer_gene_census ''' )
    cancer_gene_census_tuple = cur4.fetchall()
    df_cgc = pd.DataFrame(cancer_gene_census_tuple, columns=['Gene Symbol'])
    df_cgc =  df_cgc.rename(columns={'Gene Symbol': 'GeneName'})
    #df_cgc_all_dict = df_cgc.to_dict(orient='records')

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

    cur6 = mysql.get_db().cursor()
    cur6.execute(''' SELECT `Ensembl Transcript ID`, `Associated Transcript Name` FROM mart_export ''')
    biomart_tuple = cur6.fetchall()
    df_biomart = pd.DataFrame(biomart_tuple, columns=['Ensembl Transcript ID','Associated Transcript Name'])
    df_biomart = df_biomart.drop_duplicates()
    df_biomart =  df_biomart.rename(columns={'Ensembl Transcript ID': 'DomCancerTrans','Associated Transcript Name':'Transcript_Name'})

    result_df = pd.merge(df, df_biomart, how='inner', on='DomCancerTrans')
    result = result_df.to_dict(orient='records')


    statistic_table = df[['GeneName1', 'GeneName2', 'NumberOfStringInt', 'MissedInteractions', 'Domain1', 'Domain2', 'StringDensityRank1', 'Region1', 'DomCancerTrans']].drop_duplicates()
    statistics_table_dict = statistic_table.to_dict(orient='records')

    data_dict = df.to_dict(orient='records')


    ### DRAW PIE CHART

    data = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'},{'type':'domain'}]])

    data.add_trace(go.Pie(values=[statistic_table.iloc[0,6]*100,(100-(statistic_table.iloc[0,6]*100))], title="STRING Density Score"), 1, 1)
    data.update_traces(selector=dict(type='pie'),
                              marker=dict(colors=['lightgreen', 'darkorange'],
                              line=dict(color='#000000', width=4)))

    data.add_trace(go.Pie(labels=df.Tissue, title="Cancer Types"), 1, 2)
    data.update_traces(selector=dict(type='pie'),
                              marker=dict(line=dict(color='#000000', width=4)))

    data.update_layout(title=" Gene Name: {} ".format(genename), showlegend=False, title_font_size=24)
    graphJSON2 = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    #return render_template('network_lochen.html', genename=genename, enstid=enstid, partner_genenames=partner_genenames, data=data_dict, data_statistics = statistics_table_dict, df_cgc_list=df_cgc_list, result=result, cgc=df_cgc_dict, df_cgc_all_dict=df_cgc_all_dict, graphJSON2=graphJSON2)
    return render_template('GeneBased.html', genename=genename, data=data_dict, data_statistics = statistics_table_dict, result=result, graphJSON2=graphJSON2)

### End of Gene Based


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

@app.route("/download", methods=['GET', 'POST'])
def download():

    return render_template("download.html")


@app.route("/Sample_Based", methods=['GET', 'POST'])
def Sample_Based():

    CancerSampleId = request.args.get('sampleid')
    genename = request.args.get('gene')
    tissue = request.args.get('tissue')
   ## normal_trans_id, cancer_trans_id = encoded_img_data(CancerSampleId, genename)

    cur = mysql.get_db().cursor()
    sql = ''' SELECT Tissue, ENSG, GeneName1, CancerSampleId, DomCancerTrans, GTExMDIs FROM interactiondisruptionindominanttranscripts WHERE CancerSampleId = %s '''
    adr = (CancerSampleId, )
    cur.execute(sql, adr)
    isonet_tuple = cur.fetchall()

    df = pd.DataFrame(isonet_tuple, columns=['Tissue', 'ENSG', 'GeneName1', 'CancerSampleId', 'DomCancerTrans', 'GTExMDIs'])
    df =  df.rename(columns={'ENSG': 'ENSGid'})

    dff = df[df.GeneName1 == genename]
    cancer_trans_id = list(dff['DomCancerTrans'])[0] ## DomCancerTrans id
    normal_trans_ids = dff['GTExMDIs'].str.split(":", n = 1, expand = True)
    normal_trans_id = normal_trans_ids.iloc[0,0] # Normal Transcript Id


    return render_template("Sample_Based.html", CancerSampleId=CancerSampleId, tissue=tissue, genename=genename, normal_trans_id=normal_trans_id, cancer_trans_id=cancer_trans_id)

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

    gtex = '_gtex'
    my_tissue = tissuetype+gtex
    df_gtex = pd.read_sql(''' SELECT * FROM ''' + my_tissue,  mysql.get_db())

    pcawg = '_pcawg'
    my_tissue2 = tissuetype + pcawg
    df_pcawg = pd.read_sql(''' SELECT Feature, ''' +  '''`''' + CancerSampleId + '''`'''  + ''' FROM ''' + my_tissue2,  mysql.get_db())

    data = go.Figure()

    dff = df[df.GeneName1 == genename]
    col_name = dff.iloc[0,3] ## sample id
    cancer_trans_id = dff.iloc[0,4] ## DomCancerTrans id
    normal_trans_ids = dff['GTExMDIs'].str.split(":", n = 1, expand = True)
    normal_trans_id = normal_trans_ids.iloc[0,0] # Normal Transcript Id

    ## Transcript Counts of Transcripts in pcawg data
    count_data_cancer_trans = df_pcawg[df_pcawg.iloc[:,0].str.split(".", n = 1, expand = True)[0] == cancer_trans_id][col_name]
    count_data_normal_trans = df_pcawg[df_pcawg.iloc[:,0].str.split(".", n = 1, expand = True)[0] == normal_trans_id][col_name]

    x=[cancer_trans_id, normal_trans_id]
    y=[float(count_data_cancer_trans),float(count_data_normal_trans)]

    ## Transcript Counts of Transcripts in gtex data
    new_data_gtex = df_gtex[df_gtex.iloc[:,0].str.split(".", n = 1, expand = True)[0] == normal_trans_id]
    new_data_cancer = df_gtex[df_gtex.iloc[:,0].str.split(".", n = 1, expand = True)[0] == cancer_trans_id]

    count_data_cancer_trans_median = np.median(new_data_cancer.iloc[0,1:len(new_data_cancer.columns)])
    count_data_normal_trans_median = np.median(new_data_gtex.iloc[0,1:len(new_data_gtex.columns)])

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

    data.update_layout(title=genename, yaxis_title="TPM Count",title_font_size=24)


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