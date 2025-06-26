import db_utils
import argparse
import pandas
import numpy as np
from scipy.optimize import leastsq

'''
script to upload titration values and fits to luthy_data tables titration_values and titration_fit

3 modes of use:
- always provide -id project_id (complete name), KL_id for NL_plasmid_id and for mCit_plasmid_id (KL_...)
  -> this uploads all data related to both KL_ids within the project_id
- provide path to a tsv file containing the project_id, NL_plasmid_id and mCit_plasmid_id for all the datasets to be uploaded

example in terminal:
python titration_analysis_upload_script_v2.py -id Lu177r01 KL_177 KL_178
python titration_analysis_upload_script_v2.py -f /Users/luck/data_to_upload.txt

'''

def get_saturation_curve_data(project_id,NL_id,mCit_id,connect,FL_id):
    titration_query = f"""select a.project_id,a.NL_plasmid,a.mCit_plasmid,a.NL_property,a.mCit_property,
                c.measurement/d.measurement BRET,a.NL_plasmid_id,a.mCit_plasmid_id,a.plate_id,a.repl_id,b.avg_FL,
                e.measurement totLu, a.well_id
                from luthy_data.plate_layout a, luthy_data.FL_avg_points b, luthy_data.LU_raw c,
                    luthy_data.LU_raw d, luthy_data.LU_raw e
                where a.project_id='{project_id}' and a.include=1 and a.project_id=b.project_id and
                    a.project_id=c.project_id and a.project_id=d.project_id and a.project_id=e.project_id and
                    a.plate_id=b.plate_id and a.plate_id=c.plate_id and a.plate_id=d.plate_id and
                    a.plate_id=e.plate_id and a.well_id=b.well_id and a.well_id=c.well_id and
                    a.well_id=d.well_id and a.well_id=e.well_id and a.NL_plasmid_id='{NL_id}' and
                    a.mCit_plasmid_id='{mCit_id}' and c.measurement_id='accLu01' and d.measurement_id='donLu01'
                    and e.measurement_id='totLu01' and b.measurement_id='{FL_id}'"""
    bleedthrough_query = f"""select c.measurement/d.measurement bleedthrough,a.plate_id
                from luthy_data.plate_layout a, luthy_data.LU_raw c, luthy_data.LU_raw d
                where a.project_id='{project_id}' and a.include=1 and a.project_id=c.project_id and
                    a.project_id=d.project_id and a.plate_id=c.plate_id and a.plate_id=d.plate_id and
                    a.well_id=c.well_id and a.well_id=d.well_id and a.NL_plasmid_id='KL_03' and
                    a.mCit_plasmid_id='empty' and c.measurement_id='accLu01' and d.measurement_id='donLu01'"""
    bkg_query = f"""select b.avg_FL bkg_FL,c.measurement bkg_LU,a.plate_id
                    from luthy_data.plate_layout a, luthy_data.FL_avg_points b, luthy_data.LU_raw c
                    where a.project_id='{project_id}' and a.include=1 and a.project_id=b.project_id and
                        a.project_id=c.project_id and a.plate_id=b.plate_id and a.plate_id=c.plate_id and
                        a.well_id=b.well_id and a.well_id=c.well_id and a.NL_plasmid_id='KL_01' and
                        a.mCit_plasmid_id='empty' and b.measurement_id='{FL_id}' and c.measurement_id='totLu01'"""

    BRET_df = pandas.read_sql(titration_query,connect)
    #print(BRET_df)
    BT_df = pandas.read_sql(bleedthrough_query,connect)
#    print(BT_df)
    bkg_df = pandas.read_sql(bkg_query,connect)

    if BT_df.shape[0] == 0:
        BRET_df['cBRET'] = BRET_df['BRET'] - 0.297
        print('Warning', project_id, NL_id, mCit_id, 'have no bleedthrough control, taking default bkg value to correct.')
    else:
        BRET_df['cBRET'] = BRET_df['BRET'] - np.mean(list(BT_df['bleedthrough']))

    if bkg_df.shape[0] == 0:
        BRET_df['avg_FL'] = BRET_df['avg_FL'] - 299
        BRET_df['totLu'] = BRET_df['totLu'] - 7279
        print('Warning', project_id, NL_id, mCit_id, 'have no lumi and fluo bkg control, taking default bkg value to correct.')
    else:
        BRET_df['avg_FL'] = BRET_df['avg_FL'] - np.mean(list(bkg_df['bkg_FL']))
        BRET_df['totLu'] = BRET_df['totLu'] - np.mean(list(bkg_df['bkg_LU']))

    # set all negative avg_FL values to 0
    BRET_df.loc[BRET_df['avg_FL']<0,'avg_FL'] = 0

    # delete all rows with totLu < 15000
    BRET_df.drop(BRET_df.loc[BRET_df['totLu']<15000].index,inplace=True)

    # compute expression ratio for all remaining entries
    BRET_df['expr_ratio'] = BRET_df['avg_FL']/BRET_df['totLu']

    # set all negative cBRET ratios to 0
    BRET_df.loc[BRET_df['cBRET']<0,'cBRET'] = 0
#    print(BRET_df)

    return BRET_df


#connect to mysql db
connection = db_utils.get_connection()


#----------get input------------------------#

#initialise flags as false
parser = argparse.ArgumentParser(prog='titration_data_upload', description='uploading BRET titration data')
parser.add_argument('-id', '--id', type=str, nargs=3, help='List project_id NL_plasmid_id mCit_plasmid_id')
parser.add_argument('-f', '--input_file', type=str, nargs=1, help='Path to file containing data to be uploaded.')
args = parser.parse_args()
id = args.id
infile = args.input_file

if not id and not infile:
    raise ValueError('Error, provide information on the ids or the input file')
elif id:
    # id[0] -> project_id, id[1] -> NL_plasmid_id, id[2] -> mCit_plasmid_id
    input_specs = [id]
else:
    with open(infile[0]) as file:
        input_specs = [line[:-1].rstrip().split('\t') for line in file]

FL_id = 'FL01'

### ------  ---- ##

for dataset in input_specs:

    print(dataset)

    # fetch raw data from luthy tables
    BRET_df = get_saturation_curve_data(dataset[0],dataset[1],dataset[2],connection,FL_id)

    ###--- perform fitting ---###

    # initialize a dataframe to store the results from the fit
    BRET50_df = pandas.DataFrame({'project_id':[],'NL_plasmid_id':[],'mCit_plasmid_id':[],'repl_id':[],\
                                  'bret50':[],'bret50_err':[],'bretmax':[],'bretmax_err':[]})

    # initialize a dataframe to store all the datapoints used for making fits which then should later be uploaded to the MySQL DB
    titration_values = pandas.DataFrame({'project_id':[],'NL_plasmid_id':[],'mCit_plasmid_id':[],'repl_id':[],\
                                  'NL_property':[],'mCit_property':[],'totLu':[],'avg_FL':[],'cBRET':[]})

    repl_ids = sorted(list(BRET_df['repl_id'].unique()))

    for repl_id in repl_ids:
        sub_df = BRET_df.loc[BRET_df['repl_id'] == repl_id,['expr_ratio','cBRET','avg_FL']].reset_index()

        # only do a fit for this data if there are at least 5 datapoints and if for at least 3 of them the avg_FL values are >= 300
        if sub_df.shape[0] >= 5 and sub_df[sub_df['avg_FL']>=300].shape[0] >= 3:

            # define function to fit
            func = lambda par, adr: par[0] * adr / (par[1] + adr)
            ErrorFunc = lambda par, adr, bret: func(par, adr) - bret
            adr_list = np.array(sub_df['expr_ratio'])
            bret_list = np.array(sub_df['cBRET'])

            #tplInitial contains the "first guess" of the parameters
            Initial1 = (bret_list[len(bret_list)-1], 0.01)

            final, pcov, lsdict, lsmesg, success = leastsq(ErrorFunc, Initial1, args=(adr_list,bret_list), full_output=True)
#            print(final, pcov, lsdict, lsmesg, success)
            #final, pcov, lsdict, lsmesg, success = leastsq(ErrorFunc, Initial1, args=(adr_list, bret_list), full_output=True, ftol=0.15)
            bret_max = final[0]
            bret50 = final[1]

            # Compute std err for bret50 (according to statistical estimate of error for an estimated parameter)
            # Calculate residual variance
            # Display results
            error_values = ErrorFunc(final, adr_list, bret_list)
    #        if error_values is not None:
            s_sq = (error_values**2).sum() / (len(bret_list) - len(final))
            # else:
            #     s_sq = None

            # Multiply residual variance by fractional covariance matrix supplied by leastsq() <- this is the variance-covariance matrix
            if pcov is not None:
                pcov = s_sq * pcov
                # Extract standard error: square root of values along the diagonal of variance-covariance matrix
                sterr_bret50 = (np.absolute(pcov[1][1]))**(1/2)
                sterr_bretmax = (np.absolute(pcov[0][0]))**(1/2)
            else:
                sterr_bret50 = sterr_bretmax = 0.00

            #store the fit results
            new_df = pandas.DataFrame({'project_id':dataset[0],'NL_plasmid_id':dataset[1],'mCit_plasmid_id':dataset[2],
                                          'repl_id':repl_id,'bret50':bret50,'bret50_err':sterr_bret50,'bretmax':bret_max,
                                          'bretmax_err':sterr_bretmax},index=[0])
            BRET50_df = pandas.concat([BRET50_df,new_df])

            # store the datapoints used for the fit
            titration_values = pandas.concat([titration_values,BRET_df.loc[BRET_df['repl_id'] == repl_id,['project_id','NL_plasmid_id','mCit_plasmid_id','repl_id',\
                                          'NL_property','mCit_property','totLu','avg_FL','cBRET']]])

        else:
            print('Too few datapoints for fitting', dataset, repl_id)


    ### --- prepare and upload data to mysql --- ###


    db_table = 'luthy_data.titration_values'
    with connection.cursor() as cursor:
        query_titration_values = f"""INSERT INTO {db_table}
            (project_id, NL_plasmid_id, mCit_plasmid_id, repl_id, NL_property, mCit_property, lumi, fluo, cBRET)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE NL_property = VALUES(NL_property),mCit_property = VALUES(mCit_property),lumi = VALUES(lumi),
                fluo = VALUES(fluo), cBRET = VALUES(cBRET)"""

        for index, row in titration_values.iterrows():
            cursor.execute(query_titration_values, (row['project_id'], row['NL_plasmid_id'], row['mCit_plasmid_id'],
                                                    row['repl_id'], row['NL_property'], row['mCit_property'],
                                                    row['totLu'], row['avg_FL'], row['cBRET']))

        connection.commit()



    db_table = 'luthy_data.titration_fit'
    with connection.cursor() as cursor:
        query_titration_fit = f"""INSERT INTO {db_table}
            (project_id, NL_plasmid_id, mCit_plasmid_id, repl_id, bret50, bret50_err, bretmax, bretmax_err)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE bret50 = VALUES(bret50),bret50_err = VALUES(bret50_err),bretmax = VALUES(bretmax),bretmax_err = VALUES(bretmax_err)"""

        for index, row in BRET50_df.iterrows():
            cursor.execute(query_titration_fit, (row['project_id'], row['NL_plasmid_id'], row['mCit_plasmid_id'],
                                                    row['repl_id'], row['bret50'], row['bret50_err'],
                                                    row['bretmax'], row['bretmax_err']))

        connection.commit()

    print(dataset, 'completed')
