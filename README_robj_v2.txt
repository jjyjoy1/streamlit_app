This script used to read R expressionset object, and enable to extraact the following information.

dataset_name (string): The name of the dataset extracted from the file name without extension
X (pandas DataFrame): Expression data with samples as rows and features as columns
y (numpy array or Surv object):

If analysis is 'surv': A survival analysis object
Otherwise: Class labels as integers


groups (numpy array or None): Group identifiers if 'Group' column exists in sample metadata
group_weights (numpy array or None): Weights for groups if 'GroupWeight' column exists
sample_weights (numpy array or None): Sample weights calculated based on group counts
sample_meta (pandas DataFrame): Sample metadata from the ExpressionSet
feature_meta (pandas DataFrame): Feature metadata, including any new features added

So the scientist can track individual data shape and realted information

Upload any RDS file containing an ExpressionSet object
View the expression data matrix
Explore sample metadata
Examine feature metadata
Visualize class distributions or survival data
