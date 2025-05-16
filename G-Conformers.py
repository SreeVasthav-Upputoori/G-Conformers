import csv
import random
import multiprocessing
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, SDWriter
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField

def evaluate_energy(mol, conf_id):
    """Calculate the energy of a given conformer using UFF (Universal Force Field)."""
    ff = UFFGetMoleculeForceField(mol, confId=conf_id)
    return ff.CalcEnergy()

def process_single_smiles(smiles_population):
    """Process a single SMILES string into an optimized conformer."""
    smiles, population_size, generations, mutation_rate = smiles_population

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, "Invalid_SMILES", None, None

    mol = Chem.AddHs(mol)
    population = []

    for _ in range(population_size):
        conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        if conf_id == -1:
            continue
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
        population.append((conf_id, evaluate_energy(mol, conf_id)))

    if not population:
        return smiles, "Failed_to_generate_conformers", None, None

    population.sort(key=lambda x: x[1])

    for _ in range(generations):
        new_population = population[:population_size // 4]
        while len(new_population) < population_size:
            parents = random.sample(population[:population_size // 2], 2)
            child_id = parents[0][0]
            energy = evaluate_energy(mol, child_id)
            new_population.append((child_id, energy))
        population = sorted(new_population, key=lambda x: x[1])

    best_conf_id, best_energy = population[0]
    mol.SetProp("_Name", smiles)
    return smiles, best_energy, best_conf_id, mol

def process_smiles_parallel(smiles_list, population_size=20, generations=10, mutation_rate=0.1, num_workers=4):
    pool_inputs = [(smiles, population_size, generations, mutation_rate) for smiles in smiles_list]
    results = []

    with multiprocessing.Pool(num_workers) as pool:
        with tqdm(total=len(pool_inputs), desc="Processing", unit="mol") as pbar:
            for res in pool.imap(process_single_smiles, pool_inputs):
                results.append(res)
                pbar.update(1)

    return results

if __name__ == "__main__":
    input_file = "input_smiles.csv"
    output_file = "output_conformers.csv"
    output_sdf = "output_conformers.sdf"

    smiles_list = []
    with open(input_file, "r") as infile:
        reader = csv.reader(infile)
        header = next(reader, None)
        for row in reader:
            if row and len(row) > 0 and row[0].strip():
                smiles_list.append(row[0].strip())

    if not smiles_list:
        print("❌ Error: No valid SMILES found in input file!")
        exit()

    print(f"✅ Loaded {len(smiles_list)} SMILES molecules.")

    num_cores = min(multiprocessing.cpu_count(), len(smiles_list))
    results = process_smiles_parallel(smiles_list, num_workers=num_cores)

    valid_mols = []
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["SMILES", "Optimized Energy", "Best Conformer ID"])
        for smiles, energy, conf_id, mol in results:
            writer.writerow([smiles, energy, conf_id])
            if mol is not None and conf_id is not None:
                valid_mols.append((mol, conf_id))

    writer = SDWriter(output_sdf)
    for mol, conf_id in valid_mols:
        writer.write(mol, confId=conf_id)
    writer.close()

    print(f"\n✅ Results saved to {output_file}")
    print(f"✅ Conformers saved to {output_sdf}")
