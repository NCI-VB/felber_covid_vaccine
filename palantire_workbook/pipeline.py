# Functions defined here will be available to call in
# the code for any table.
from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.5b78a621-1c77-4be9-bed7-e2c6cc905a96")
)
from pyspark.sql.types import *
def All_Genes():
    schema = StructType([StructField("Gene", StringType(), True)])
    return spark.createDataFrame([["IL-6"],["IFN-gamma"],["IL-8"],["TNF-alpha"],["IL-12IL-23p40"],["IL-15"],["IL-16"],["IL-7"],["VEGF-A"],["CRP"],["ICAM-1"],["SAA"],["IL-17C"],["IL-17B"],["IL-1Ra"],["IL-3"],["TSLP"],["Eotaxin"],["Eotaxin-3"],["IP-10"],["MCP-1"],["MCP-4"],["MDC"],["MIP-1alpha"],["MIP-1beta"],["TARC"],["VCAM-1"]], schema=schema)

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.3097ec7a-49e0-409f-816e-7946edc18aed")
)
from pyspark.sql.types import *
def partial_patients():
    schema = StructType([StructField("patient_id", StringType(), True)])
    return spark.createDataFrame([["V150"],["V151"],["V153"],["V155"]], schema=schema)

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.4ba63091-d7c1-4e82-b2b9-6fc18808feab")
)
from pyspark.sql.types import *
def selected_analytes():
    schema = StructType([StructField("Gene", StringType(), True)])
    return spark.createDataFrame([["CRP"],["Eotaxin"],["IFN-gamma"],["IL-12IL-23p40"],["IL-15"],["IL-16"],["IL-1Ra"],["IL-3"],["IL-6"],["IL-7"],["IL-8"],["IP-10"],["MIP-1alpha"],["MIP-1beta"],["SAA"],["TNF-alpha"],["MCP-1"],["MDC"],["VEGF-A"]], schema=schema)

