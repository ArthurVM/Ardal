from ardal import Ardal

data = ["/home/arthur/BioInf/Ardal_MAIN/Ardal/tests/data/sim_matrix.npz", "/home/arthur/BioInf/Ardal_MAIN/Ardal/tests/data/sim_headers.json"]
# data = ["/home/arthur/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps/cog2M.bin", "/home/arthur/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps/cog2M.json"]
ard = Ardal(data, roaring="true", allele_id_format="{ref}.{start}.{alt}", quiet_init=False, verbosity="debug")

ard.distance.pairwise(backend="auto", metric="jaccard", as_dataframe=False, threads=12)
ard.distance.pairwise(backend="auto", metric="hamming", as_dataframe=False, threads=12)
ard.distance.pairwise(backend="auto", metric="inner_product", as_dataframe=False, threads=12)
ard.distance.neighbourhood("GUID1", n=10)
ard.distance.snv_neighbourhood("GUID1", n=10)