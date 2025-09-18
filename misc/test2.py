from ardal import Ardal

# data = ["/home/arthur/BioInf/Ardal_MAIN/Ardal/tests/data/sim_matrix.npz", "/home/arthur/BioInf/Ardal_MAIN/Ardal/tests/data/sim_headers.json"]
data = ["/home/arthur/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps/cog2M.bin", "/home/arthur/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps/cog2M.json"]
ard = Ardal(data, roaring=False, allele_id_format="{chr}.{start}.{ref}.{alt}", quiet_init=False, verbosity="debug")

ard.distance.pairwise(backend="auto", metric="jaccard", as_dataframe=False, threads=12)