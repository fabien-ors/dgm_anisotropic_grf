#%%
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import cholesky, solve_triangular
from sklearn.gaussian_process.kernels import Matern
import click
import h5py
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import os

def extract_samples(indices_to_extract,path_data,data_name):
	id_sort=np.argsort(indices_to_extract)
	id_sorted=np.array(indices_to_extract)[id_sort]
	with h5py.File(path_data, 'r') as f:
		dataset = f[data_name][id_sorted]
		if len(dataset.shape)>3:
			return dataset[:,0,:,:][id_sort]
		else:
			return dataset[id_sort]

def get_nsamples(path_data,data_name):
	with h5py.File(path_data, 'r') as f:
		return f[data_name].shape[0]

#%%

def extract_ref_data(ref_data_path):
	with h5py.File(ref_data_path, 'r') as f:
		dat_all=f["data"][0,0]
		mask=f["data"][0,1]
		dat_obs=f["data"][0,2]
	id_obs=np.where(mask.flatten()==1)[0]
	id_pred=np.where(mask.flatten()==0)[0]
	return dat_all, mask, dat_obs, id_obs, id_pred


#%%

## Create Matern kernel, ν=2, σ²=1
Mtn=Matern(length_scale=1,nu=2)
def matern_nu2(X):
    return Mtn.__call__(2*X) ## Factor 2 in accordance with parametrization used in paper

## Function to "cancel out" anisotropy
def anisotropic_change(X, length_scales, theta):
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    L = np.diag(1.0 / np.array(length_scales))
    A =  (L @ R.T) @ X.T
    return A.T

# Matern log-likelihood
def log_likelihood_matern(params,y,coords):
    n_points=len(y)
    a, rho, theta = params
    ℓ1, ℓ2 = a, a*rho
    if ℓ1 <= 0 or ℓ2 <= 0 or not (0 <= theta <= np.pi):
        return -np.inf
    try:
        D = anisotropic_change(coords, [ℓ1, ℓ2], theta)
        K = matern_nu2(D) 
        L = cholesky(K, lower=True)
        alpha = solve_triangular(L.T, solve_triangular(L, y, lower=True))
        ll = -0.5 * y @ alpha - np.sum(np.log(np.diag(L))) - n_points / 2 * np.log(2 * np.pi)
        return ll
    except np.linalg.LinAlgError:
        return -np.inf

#%%

## Function to compute MLE
def compute_MLE(y,X,initial_guesses,bounds):

	def min_f(params):
		return -log_likelihood_matern(params,y,X)

	res=[]
	for i in range(len(initial_guesses)):
		initial_guess = initial_guesses[i]  
		result = minimize(min_f, initial_guess, bounds=bounds,method="L-BFGS-B")
		if result.success:
			est_a, est_rho, est_theta = result.x
			res.append([est_a, est_rho, est_theta,-result.fun])
	res=np.array(res)
	res=res[res[:,3].argmax()]
	return res


#%%

# Worker function
def process_sample(i, dataset, Xall, Ngd, nobs, initial_guesses, bounds,temp_res_folder,id_sample,id_pred):
	print("Sample", i + 1, "/", len(dataset))
	id_obs = np.random.choice(id_pred, nobs, replace=False)
	X = Xall[id_obs]
	y = dataset[i].flatten()
	y = y[id_obs]
	res=compute_MLE(y, X, initial_guesses, bounds)
	if temp_res_folder is not None:
		np.savetxt(temp_res_folder+"res_"+str(id_sample)+".csv",res,delimiter=",",header="a,rho,theta,llmax")
	return res

#%%

@click.command()
@click.option(
    "--data_path",
    help="name of samples folder",
    metavar="STR",
    type=str,
    required=True
)
@click.option(
    "--ref_data_path",
    help="name of data folder with mask",
    metavar="STR",
    type=str,
    default=None,
    required=False
)
@click.option(
    "--seed",
    help="seed",
    metavar="INT",
    type=int,
    default=42,
    required=False
)
@click.option(
    "--nsamples",
    help="max number of samples",
    metavar="INT",
    type=int,
    default=100,
    required=False
)
@click.option(
    "--id_s0",
    help="first index of samples",
    metavar="INT",
    type=int,
    default=0,
    required=False
)
@click.option(
    "--nobs",
    help="number of observations",
    metavar="INT",
    type=int,
    default=1000,
    required=False
)
@click.option(
    "--data_name",
    help="name of samples dataset in h5 file",
    metavar="STR",
    type=str,
    default="images",
    required=True
)
@click.option(
    "--file_name",
    help="export name",
    metavar="STR",
    type=str,
    default="export.csv",
    required=False
)
@click.option(
    "--nworkers",
    help="number of workers",
    metavar="INT",
    type=int,
    default=5,
    required=False
)
@click.option(
    "--temp_res_folder",
    help="folder to save intermediate results",
    metavar="STR",
    type=str,
    default=None,
    required=False
)
def cmdline(
    data_path: str,
    data_name:str,
    seed: int,
    nsamples: int,
    nobs: int,
    file_name: str,
    nworkers: int,
    temp_res_folder: str,
    id_s0: int,
    ref_data_path: str,
    **opts,
):
	
	## Bounds for paramters
	bounds = [(0.05, 0.3), (0.1, 1), (0, np.pi)]
	## Initial guess for optimization
	initial_guesses=[[0.15,0.5,np.pi/2], [0.075,0.2,np.pi/2]]

	## Load data
	ntot_im=get_nsamples(data_path,data_name)
	print("Total number of images:",ntot_im)

	## Create grid coordinates
	Ngd=256
	Xgd=np.meshgrid(np.linspace(0,1,Ngd),np.linspace(0,1,Ngd))
	Xall=np.c_[(Xgd[0][np.flip(range(Ngd))].flatten(),Xgd[1][np.flip(range(Ngd))].flatten())]

	## Extract mask
	if ref_data_path is None:
		ref_data_path=data_path
	dat_all, mask, dat_obs, id_obs, id_pred = extract_ref_data(ref_data_path)

	## Select samples
	id_s=np.array(range(id_s0,id_s0+nsamples))
	dataset=extract_samples(id_s,data_path,data_name)

	## Tenporary result folder
	if temp_res_folder is None:
		temp_res_folder="temp_"+(data_path.split("/")[-1].split(".")[0])+"/"
	if not os.path.exists(temp_res_folder):
		os.makedirs(temp_res_folder)

	# Run
	results_list=list()
	for i in range(len(id_s)):
		if np.any(np.isnan(dataset[i])):
			continue
		res=process_sample(i, dataset, Xall, Ngd, nobs, initial_guesses, bounds,temp_res_folder,id_s[i],id_pred)
		results_list.append(res)
	# Save results
	results = np.array(results_list)
	np.savetxt(file_name,results,delimiter=",",header="a,rho,theta,llmax")

if __name__ == "__main__":
    cmdline()


