import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC,SVC
from sklearn.decomposition import PCA
from tqdm import tqdm
import pdb

def main():

    #Define output 
    output_dir = "./analysis/linearSVM_out"

    #data paths
    counts_path = "./analysis/normalized_counts/Allcounts_batchnorm_vst.csv"
    metadata_path = "./data/Example_data/Example_metadata_Ranalysis.tsv"
    kmeans_CBC_path = "./analysis/kmeans/kmeans_groups_cbc.csv" 
    kmeans_all_path = "./analysis/kmeans/kmeans_groups_all.csv" 

    validation_counts_path = "./data/validation_data/ALL_counts_batchnorm_vst.csv"
    validation_metadata_path = "./data/validation_data/ALL_metadata.csv"

    #read in the data 
    count_data = pd.read_csv(counts_path, index_col=0, sep = ",")
    metadata = pd.read_csv(metadata_path, sep = "\t")
    kmeans_data_cbc = pd.read_csv(kmeans_CBC_path)
    kmeans_data_all = pd.read_csv(kmeans_all_path)

    validation_counts = pd.read_csv(validation_counts_path, index_col=0, sep = ",")
    validation_metadata = pd.read_csv(validation_metadata_path, sep = ",")

    #Feels bad to write this...need to fix #'s in SAW names later
    validation_counts.columns=validation_counts.columns.str.replace('.','#') 

    #Merge kmeans with metdata
    metadata = metadata.merge(kmeans_data_cbc, on = "Sample_ID", how = "left")
    metadata = metadata.merge(kmeans_data_all, on = "Sample_ID", how = "left")

    #Train on CBC data and test on samples from other groups 
    #categories = ["IDSA_SCORE_1to4", "Ulcer_duration_cat", "cluster_res_cbc"]
    categories = ["cluster_res_cbc", "IDSA_SCORE_1to4", "Ulcer_duration_cat"]

    for cat in categories:

        cleaned_data = clean_data(count_data, metadata, "Sample_ID", cat)
        print(cleaned_data[0].shape)
        print(cleaned_data[1].shape)

        # plot_nFeatVsAcc(cleaned_data, total_features= 100, cat = cat, count_data = count_data, \
        #                 metadata = metadata, save=True)

        model_featureselect = Linear_classifier(cleaned_data[0], cleaned_data[1], classifier_term = cat, \
                                                subset = "Source == 'CBC'", max_features = 20)
        genes = model_featureselect.genes
        good_genes = model_featureselect.good_genes

        print(f"For the predictor {cat}, the good genes are: {good_genes}")
      
        # if cat == "cluster_res_cbc": 
        #     cat = "cluster_res_all"

        ######## Model Training ######### 
        train_data = clean_data(count_data, metadata, "Sample_ID", cat, subset = "Source == 'CBC'")
        test_data =  clean_data(count_data, metadata, "Sample_ID", cat, subset = "Source != 'CBC'")

        #Mask non-useful genes from train/test data
        mask = np.array([1 if x in good_genes else 0 for x in genes]) == 1
        train_data[0] = np.array(train_data[0])[mask, :]  
        test_data[0] = np.array(test_data[0])[mask, :] 
        model_testtrain = SVC(C=0.8, max_iter=10000, probability = True)
        #model_testtrain = LinearSVC(C=1, penalty='l2', dual=False, max_iter=10000)
        model_testtrain.fit(train_data[0].T, np.array(train_data[1][cat], dtype = np.int))

        # #test
        # accuracy = model_testtrain.score(test_data[0].T, np.array(test_data[1][cat], dtype = np.int))
        # print(test_data[0].shape)

        # train_accuracy = model_testtrain.score(train_data[0].T, np.array(train_data[1][cat], dtype = np.int))

        # print(f"The test accuracy for {cat} was: {accuracy}")
        # print(f"The train accuracy for {cat} was: {train_accuracy}")

        #Test classification for validation data for cat = cluster_res_cbc:
        if cat == "cluster_res_cbc":

            validation_counts = validation_counts.loc[good_genes]

            validation_data = clean_data(validation_counts, validation_metadata, "Sample_ID", "Specific_ID", \
                subset = "Source != 'CBC'")

            print([data.shape for data in validation_data])

            print(model_testtrain.predict(validation_data[0].T))

            prediction = model_testtrain.predict(validation_data[0].T)
            prediction_probability = model_testtrain.predict_proba(validation_data[0].T) 

            results = pd.DataFrame({'Sample_ID':validation_data[1]["Sample_ID"], \
                                    'prediction': prediction.tolist(),\
                                    'prediction_proba:' : prediction_probability.tolist()})

            results.to_csv("./data/validation_data/predictions.csv", index = False)
        
        #Export Results
        for i in range(0, len(model_featureselect.lr.classes_)):
            cat_name = model_testtrain.classes_[i] #So that the names are correct           
            out_goodgenes = pd.DataFrame(good_genes, columns = ["Predictor Genes"])
            out_goodgenes_filename = f"{output_dir}/GoodGenes_{cat}.csv"
            out_goodgenes.to_csv(out_goodgenes_filename, index = False)

        #If using linearSVC, can remove comments and get coefficients/plots 
        for i in range(0,len(model_featureselect.lr.coef_)):
            cat_level = model_testtrain.classes_[i]
            model_coefs = model_featureselect.lr.coef_[i]
            pred_genes = model_featureselect.good_genes

            # Mask unused coefs
            coef_mask = np.array([1 if x in good_genes else 0 for x in genes]) == 1
            masked_coefs = model_featureselect.lr.coef_[i,coef_mask]
                        
            title = f"{cat}: {cat_level}"
           
            plot_fimportances(masked_coefs, good_genes, cat=title)
            
            #Output the data:
            out_data = np.column_stack((pred_genes, masked_coefs))
            out_file = f"{output_dir}/COEF_{cat}_{cat_level}.csv"
            out_data = pd.DataFrame(out_data, columns = ["Genes", "Coefficient"])
            out_data = out_data.sort_values("Coefficient", ascending = False)
            out_data.to_csv(out_file, index = False)

def clean_data(data, metadata, sample_id_colname, classifier_term, subset = False ):
    """ This function will take a rna-seq count data set, metadata,
        and identifying col name to remove samples from metadata that aren't
        in the main data. 

        It returns a list where the first element is the metadata and the second is
        the metadata.

        It will also make sure the metadata and count data are in the proper order 

        Note: If classifier term, it returns the data as np arrays instead of pandas(actally doesnt 
        look like this is true )
    """

    #Remove samples from metadata which aren't in the data
    metadata = metadata.loc[metadata[sample_id_colname].isin(data.columns)]
    data = data.loc[:, data.columns.isin(metadata[sample_id_colname])]

    if subset: 
        metadata = metadata.query(subset)
    
    if classifier_term:
        metadata = metadata.query(f"{classifier_term} != 'NA' & {classifier_term} != 'NaN'" )
        data = data.loc[:, data.columns.isin(metadata[sample_id_colname])]

    #Fix order or data and metadata 
    data = data.reindex(sorted(data.columns), axis = 1)
    metadata = metadata.sort_values(by=["Sample_ID"])

    cols_equal = np.array_equal(data.columns, metadata["Sample_ID"])

    if cols_equal:
        return([data, metadata])

    else:
        print("ERROR: Orders of samples in metadata/data don't match")
        sys.exit()

def plot_fimportances(coef, f_names, cat):
    imp = coef
    imp,f_names = zip(*sorted(zip(imp,f_names))) #Added 0 to imp and it worked
    plt.barh(range(len(f_names)), imp, align='center')
    plt.yticks(range(len(f_names)), f_names)
    plt.title(cat)
    plt.tight_layout(pad=1)
    #plt.show()

def plot_nFeatVsAcc(cleaned_data, count_data, metadata, total_features, cat, save = False):

    print("Now Testing the accuracy for different #'s of features")

    model_accuracy = []
    feature_nr = []

    #Run the models and get the accuracy
    for n_feature in tqdm(range(1,total_features)):
        cleaned_data = clean_data(count_data, metadata, "Sample_ID", cat)
        #print(cleaned_data[0].shape)
        #rint(cleaned_data[1].shape)
        model_featureselect = Linear_classifier(cleaned_data[0], cleaned_data[1], classifier_term = cat, \
                                                subset = "Source == 'CBC'", max_features = n_feature)
        genes = model_featureselect.genes
        good_genes = model_featureselect.good_genes

        #print(f"For the predictor {cat}, the good genes are: {good_genes}")

        if cat == "cluster_res_cbc": 
            cat = "cluster_res_all"

        ######## Model Training ######### 
        train_data = clean_data(count_data, metadata, "Sample_ID", "cluster_res_cbc", subset = "Source == 'CBC'")
        test_data =  clean_data(count_data, metadata, "Sample_ID", "cluster_res_all", subset = "Source != 'CBC'")

        #Mask non-useful genes from train/test data
        mask = np.array([1 if x in good_genes else 0 for x in genes]) == 1
        train_data[0] = np.array(train_data[0])[mask, :]  
        test_data[0] = np.array(test_data[0])[mask, :] 
        #model_testtrain = SVC(C=1, max_iter=10000)
        model_testtrain = LinearSVC(C=0.1, penalty='l2', dual=False, max_iter=10000)
        model_testtrain.fit(train_data[0].T, np.array(train_data[1][cat], dtype = np.int))
        
        #test
        accuracy = model_testtrain.score(test_data[0].T, np.array(test_data[1][cat], dtype = np.int))
        train_accuracy = model_testtrain.score(train_data[0].T, np.array(train_data[1][cat], dtype = np.int))

        model_accuracy.append(accuracy)
        feature_nr.append(n_feature)

    #Make the figure

    fig, ax = plt.subplots()
    ax.plot(feature_nr, model_accuracy)
    ax.set(xlabel = 'Number of Genes in Classifier', ylabel = 'Classifier Accuracy')
    #plt.show()

    if save:
        fig.savefig(f"./analysis/Figures/Accuracy_{cat}.png", dpi = 300, bbox_inches = 'tight')

class Linear_classifier():
    """ Class for running the linear classifier
    data_path - Path to count data. Data should be a data frame of normalized counts exported from R as csv
                where the first row contains the sample_IDs and the first column contains the genes included
                in the analysis.

    metadata_path - a csv file containg the sample ids and any relevant metadata

    classifier_term - a string the metadata variable used to train the classifier against

    subset - a pandas.query() string which is used to subset the data [optional]
    """
    def __init__(self, data, metadata, classifier_term, max_features , subset=False):
        #self.data_path = data_path
        #self.metadata_path = metadata_path
        self.data = data
        self.metadata = metadata
        self.subset = subset    
        self.classifier_term = classifier_term
        self.max_features = max_features

        #run the analysis
        self.prepare_data()
        self.run_analysis()
        self.make_plots()

    def prepare_data(self):

        self.genes = self.data.index.values

        if self.subset:
            #Not sure if query is best option here, but seems to work ok
            self.metadata = self.metadata.query(self.subset)
            #Remove samples from count data if not in metadata
            sample_ids = self.metadata["Sample_ID"]           
            #Remove samples from data if not in subsetted metadata  
            self.data = self.data.loc[:, self.data.columns.isin(sample_ids)]

        #Remove from metadata if not in actual count data
        self.sample_ids = self.data.columns
        self.metadata = self.metadata[self.metadata["Sample_ID"].isin(self.sample_ids)]

    def run_analysis(self):

        #Generate PCA - still not sure why this is different from the R PCA output
        pca = PCA(n_components=2)
        self.pcacomps = pca.fit_transform(self.data.T)
        
        #Get a mask index for NaNs
        NA_mask = ~np.isnan(self.metadata[self.classifier_term])
        model_data = np.array(self.data)
        #print("Before NaN removal: ",model_data.shape)
        
        #NaN values are filtered
        model_data = model_data[:, NA_mask]
        #print("After NaN removal: ", model_data.shape)

        #Get classifier term with na.removed
        classifier_term = self.metadata[self.classifier_term][NA_mask]

        #Fit the linear classifier
        self.lr = LinearSVC(C=0.05, penalty='l1', dual=False, max_iter=10000)
        self.lr.fit(model_data.T, np.array(classifier_term, dtype=np.int))

        #Get important genes
        #print('Before Feature Selection =', model_data.shape)
        model = SelectFromModel(self.lr, prefit=True, threshold=-np.inf, max_features=self.max_features)
        #print(model)
        #print('After Feature Selection =', model.transform(model_data.T).T.shape)
        self.good_genes =  self.genes[model.get_support()]
    #    print(f"For {self.classifier_term}, these genes are good: {self.good_genes}")

    def make_plots(self):

        # plt.subplot(2, 2, 1, label='PCA')
        # plt.scatter(self.pcacomps[:, 0], self.pcacomps[:, 1], c = self.metadata[self.classifier_term])

        c = self.lr.coef_
        proj = c @ self.data

        # plt.subplot(2, 2, 2, label='01')
        # plt.scatter(proj.iloc[0, :], proj.iloc[1, :], c=self.metadata[self.classifier_term])
        # plt.subplot(2, 2, 3, label='12')
        # plt.scatter(proj.iloc[1, :], proj.iloc[2, :], c=self.metadata[self.classifier_term])
        # plt.subplot(2, 2, 4, label='02')
        # plt.scatter(proj.iloc[0, :], proj.iloc[2, :], c=self.metadata[self.classifier_term])

        #plt.show()


if __name__ == "__main__": 
    main()
