\name{scMAGIC}
\alias{scMAGIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A cell type classifier for single cell RNA sequencing data
%%  ~~function to do ... ~~
}
\description{
scMAGIC main function. User can provide reference data and query data,
and scMAGIC can assign labels to cells in query dataset based on reference.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scMAGIC(exp_sc_mat, exp_ref_mat, exp_ref_label = NULL,
        single_round = F, identify_unassigned = T,
        atlas = c(NULL, 'MCA', 'HCL'), use.RUVseq = T,
        cluster_num_pc = 50, cluster_resolution = 3,
        combine_num_cell = NULL, min_cell = 1,
        method1 = ('kendall', 'spearman', 'pearson', 'cosine', 'multinomial'),
        method2 = ('multinomial', 'kendall', 'spearman', 'pearson', 'cosine', 'randomforest'),
        corr_use_HVGene1 = 2000, corr_use_HVGene2 = 2000,
        threshold = 5, num_threads = 4, cluster_assign = F, simple.output = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp_sc_mat}{The expression matrix of query data}
  \item{exp_ref_mat}{The expression matrix of reference data. If the type_ref is 'sc-counts', 'exp_ref_mat' is a single cell counts matrix where the row names correspond to gene symbols and the column names correspond to cell barcodes; if not, 'exp_ref_mat' is a matrix where the row names correspond to gene symbols and the column names correspond to cell type labels.}
  \item{exp_ref_label}{If 'type_ref' is 'sc-counts', 'exp_ref_label' is the vector including well-annotated cell labels corresponding to the cell barcodes of the 'exp_ref_mat'; if not, 'exp_ref_label' is NULL.}
  \item{single_round}{Whether to use single round annotation strategy, by default, scMAGIC use two-round annotation strategy.}
  \item{identify_unassigned}{Whether to label some cells as "Unassigned", by default, it is TRUE. If you think reference covers all cell types in query data, choose FALSE.}
  \item{atlas}{If reference is from mouse, it is 'MCA'; if reference is from human, it is 'HCL'; if no reference is avaiable, it is NULL.}
  \item{use_RUVseq}{Whether to use 'RUVSeq' to remove the batch effect between reference and atlas, by default, choose TRUE.}
  \item{cluster_num_pc}{Number of PCs used in clustering.}
  \item{cluster_resolution}{Resolution of the clustering algorithm, if you want get more clusters, you can give a larger resolution.}
%%  \item{opt_speed}{Whether to speed calculation by merge the expression profiles of similar cells.}
  \item{combine_num_cell}{Mean combined cell number in speeding calculation.}
  \item{min_cell}{If the number of validated cells with a reference cell type isn't less than "min_cell", these cells will be added into local reference.}
  \item{method1}{The methods of similarity calculation in first-round annotation, by default, choose 'kendall'.}
  \item{method2}{The methods of similarity calculation in second-round annotation, by default, choose 'multinomial'.}
  \item{corr_use_HVGene1}{Number of genes to select as top variable genes in similarity calculation of first-round annotation}
  \item{corr_use_HVGene2}{Number of genes to select as top variable genes in similarity calculation of second-round annotation}
  \item{threshold}{If confidence score is lower than 'threshold', we consider the corresponding label incorrect.}
  \item{num_threads}{Number of CPU used in calculation.}
  \item{cluster_assign}{Whether to annotate single cells by assigning a cell type to each cluster.}
  \item{simple_output}{Whether to output intermediate results.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe including cell type labels and confidence scores.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yu Zhang
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
