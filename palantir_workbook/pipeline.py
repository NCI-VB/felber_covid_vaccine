# Functions defined here will be available to call in
# the code for any table.
from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.4ba63091-d7c1-4e82-b2b9-6fc18808feab")
)
from pyspark.sql.types import *
def selected_analytes():
    schema = StructType([StructField("Gene", StringType(), True)])
    return spark.createDataFrame([["CRP"],["Eotaxin"],["IFN-gamma"],["IL-12IL-23p40"],["IL-15"],["IL-16"],["IL-1Ra"],["IL-3"],["IL-6"],["IL-7"],["IL-8"],["IP-10"],["MIP-1alpha"],["MIP-1beta"],["SAA"],["TNF-alpha"],["MCP-1"],["MDC"],["VEGF-A"]], schema=schema)

