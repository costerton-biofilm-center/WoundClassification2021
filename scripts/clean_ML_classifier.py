import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.svm import SVC, LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import cross_val_score, StratifiedKFold, ShuffleSplit, KFold
from tqdm import tqdm
import pdb 

def main():
	# Output dir 
	output_dir = "./analysis/linearSVM_out"

	# Data Paths
	counts_path = "./analysis/counts/Batchnorm_mRNA_vst.csv"
	metadata_path = "./analysis/metadata/metadata_with_results.csv"

	# Get the data 
	count_data = pd.read_csv(counts_path, index_col=0, sep = ",")
	metadata = pd.read_csv(metadata_path, sep = ",")

	#Filter metdata, format variables, make dummies 
	metadata = metadata[["Sample_ID", "Source", "cluster_res_all", "high_bacteria_1High0Low", "PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]]
	dummies = pd.get_dummies(metadata["cluster_res_all"], prefix = "isCluster")

	#Combine the data'
	metadata = pd.concat([metadata, dummies], axis=1)
	
	#Select training categories
	categories = ["cluster_res_all", "PEDIS_IDSA_1uninfected_2mild_3mod_4severe"]

	#Feature Selection 

	for cat in categories: 
		#Make sure the data is clean & only from CBC 
		cleaned_data = clean_data(count_data, metadata, "Sample_ID", cat)

		#Fit the model 
		model = LinearSVC(C=1, penalty='l1', dual=False, max_iter=100000)
		model.fit(np.array(cleaned_data[0]).T, np.array(cleaned_data[1][cat], dtype = np.int)) 
		
		#Get the good predictor genes
		model_fselect = SelectFromModel(model, prefit=True, threshold=-np.inf, max_features=20)
		good_genes = count_data.index[model_fselect.get_support()]
		
		print(f"For {cat}, the good genes are : {good_genes}")

		plot_nFeatVsAcc(cleaned_data[0], cleaned_data[1], cat, n_features = 100, save = True)

		#Make a mask for good genes 
		mask = np.array([1 if x in good_genes else 0 for x in count_data.index]) == 1

		#Export the results
		for i in range(0, len(model.classes_)):
			cat_name = model.classes_[i] #So that the names are correct
			out_goodgenes = pd.DataFrame(good_genes, columns = ["Predictor Genes"])
			out_goodgenes_filename = f"{output_dir}/GoodGenes_{cat}.csv"
			out_goodgenes.to_csv(out_goodgenes_filename, index = False)

		#If using linearSVC, can remove comments and get coefficients/plots 
		for i in range(0,len(model.coef_)):
			cat_level = model.classes_[i]
			model_coefs = model.coef_[i]

			#Mask unused coefs
			masked_coefs = model.coef_[i,mask]
            
			#Define title and make plots
			title = f"{cat}: {cat_level}"
			plot_fimportances(masked_coefs, good_genes, cat=title)
            
			#Output the data:
			out_data = np.column_stack((good_genes, masked_coefs))
			out_file = f"{output_dir}/COEF_{cat}_{cat_level}.csv"
			out_data = pd.DataFrame(out_data, columns = ["Genes", "Coefficient"])
			out_data = out_data.sort_values("Coefficient", ascending = False)
			out_data.to_csv(out_file, index = False)

		# # Train and fit using only good genes
		# valid_data_train = clean_data(validation_counts, metadata[metadata['Source'] == 'CBC'], "Sample_ID", cat)
		# valid_data_test = clean_data(validation_counts, metadata[metadata['Source'] != 'CBC'], "Sample_ID", cat)

		# model_validation = SVC(probability=True)
		# model_validation.fit(np.array(valid_data_train[0][mask]).T, np.array(valid_data_train[1][cat], dtype = np.int)) 

		# # Predict on test data
		# test_samples = metadata['Sample_ID'][metadata['Data_set'] == 'Validation']
		# test_counts = valid_data_test[0].filter(items = test_samples)[mask]

		# prediction = model_validation.predict(np.array(test_counts).T)
		# prediction_prob = model_validation.predict_proba(np.array(test_counts).T)

		# results = pd.DataFrame({'Sample_ID':test_counts.keys(), \
		# 						'prediction': prediction.tolist(),\
		# 						'prediction_proba' : prediction_prob.tolist()})
        
		# results.to_csv(f"./analysis/validation/{cat}_predictions.csv", index = False)


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

def plot_nFeatVsAcc(counts, metadata, cat, n_features=10, save = False):

	print("Now Testing the accuracy for different #'s of features")

	model_accuracy = []
	feature_nr = []

	#Run the models and get the accuracy
	for n_feature in tqdm(range(1,n_features)):
		# Get Good Genes
		model = LinearSVC(C=1, penalty='l1', dual=False, max_iter=100000, random_state = 15815)
		model.fit(np.array(counts).T, np.array(metadata[cat], dtype = np.int)) 
		model_fselect = SelectFromModel(model, prefit=True, threshold=-np.inf, max_features=n_feature)
		good_genes = counts.index[model_fselect.get_support()]

		#Mask non-useful genes from train/test data
		mask = np.array([1 if x in good_genes else 0 for x in counts.index]) == 1

		#Fit training data 
		model = SVC(C=1, random_state = 15815)    
		
		# Split test and train 
		skf = StratifiedKFold(n_splits = 2)

		#Cross validation
		scores = cross_val_score(model, counts[mask].T, metadata[cat], cv = skf)
		
		model_accuracy.append(np.mean(scores))
		feature_nr.append(n_feature)

    #Make the figure
	fig, ax = plt.subplots()
	ax.plot(feature_nr, model_accuracy)
	ax.set(xlabel = 'Number of Genes in Classifier', ylabel = 'Classifier Accuracy')
	ax.set_xlim(0,50)
	ax.set_ylim(0.25,1.1)
	plt.yticks(np.arange(0,1.2,0.2))
	#plt.show()
	
	if save:
		fig.savefig(f"./analysis/Figures/Accuracy_{cat}.png", dpi = 300, bbox_inches = 'tight')

if __name__ == "__main__":
	main()