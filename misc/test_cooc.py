from ardal import Ardal
import time

# data = ["/home/amorris/BioInf/Ardal_MAIN/cw_ardal_comparison_WD/sim_matrix.npy", "/home/amorris/BioInf/Ardal_MAIN/cw_ardal_comparison_WD/sim_headers.json"]
data = ["/home/amorris/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps.1-5000.npz", "/home/amorris/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps.1-5000_headers.json"]
# data = ["/home/amorris/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps.1-100000.npz", "/home/amorris/BioInf/Ardal_MAIN/Ardal_WD/cog_data/cog_snps.1-100000_headers.json"]
# data = "/home/amorris/BioInf/Crypto_popgen/data/concat_crypto.csv"

ard = Ardal(data, force_roaring=True)

cooc = ard.alleleCooccurrance(threshold=0.95, use_simd=False, cache_bytes=100000000)  # 1 GB cache
print(cooc)