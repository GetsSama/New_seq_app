#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 18:05 2022

@author: SmirnygaTotoshka
"""


import argparse
import json
import os
import traceback
import multiprocessing
import time
import sys
import pandas as pd
from rdkit import Chem

alphabets = {"protein" : list("ACDEFGHIKLMNPQRSTVWY"),
             "DNA" : list("ATGC"),
             "RNA" : list("AUGC")}

def isCorrectSequence(sequence, alphabet):
    seq = sequence.strip()
    for i in range(0,len(seq)):
        if seq[i] not in alphabet:
            return False
    return True

def chargePeptide(mol):
    carboxyl = Chem.MolFromSmarts('C(=O)O')# aspartate and glutamate
    lysine = Chem.MolFromSmarts('C(CCN)C[C@@H](C(=O))')
    imydazoline = Chem.MolFromSmarts('c1cnc[nH]1') # histidine
    guanidine = Chem.MolFromSmarts('NC(N)=N') # arginine

    mol.GetAtomWithIdx(0).SetFormalCharge(1) # charge N end of peptide

    carboxyl_groups_pos = list(mol.GetSubstructMatches(carboxyl))

    if len(carboxyl_groups_pos) != 0:
        for g in carboxyl_groups_pos:
            atoms_pos = list(g)
            mol.GetAtomWithIdx(atoms_pos[2]).SetFormalCharge(-1)

    lysine_groups_pos = list(mol.GetSubstructMatches(lysine))

    if len(lysine_groups_pos) != 0:
        for g in lysine_groups_pos:
            atoms_pos = list(g)
            mol.GetAtomWithIdx(atoms_pos[3]).SetFormalCharge(1)

    imydazoline_groups_pos = list(mol.GetSubstructMatches(imydazoline))

    if len(imydazoline_groups_pos) != 0:
        for g in imydazoline_groups_pos:
            atoms_pos = list(g)
            mol.GetAtomWithIdx(atoms_pos[4]).SetFormalCharge(1)

    guanidine_groups_pos = list(mol.GetSubstructMatches(guanidine))

    if len(guanidine_groups_pos) != 0:
        for g in guanidine_groups_pos:
            atoms_pos = list(g)
            mol.GetAtomWithIdx(atoms_pos[3]).SetFormalCharge(1)

def write(proc, partition, output, sequence_column, mutation_column, isCharged, alphabet, output_filename):
    out_logs = output + "/logs"
    if not os.path.exists(out_logs):
        os.mkdir(out_logs)
    out = os.path.join(out_logs, output_filename + "_thread_" + str(proc) + ".sdf")
    log = os.path.join(out_logs, output_filename + "_thread_" + str(proc) + "_log" + ".txt")
    fail = os.path.join(out_logs, output_filename + "_thread_" + str(proc) + "_failed" + ".txt")
    try:
        with open(out, "w", encoding="utf-8") as o, open(log, "w", encoding="utf-8") as l, open(fail, "w", encoding="utf-8") as f:
            l.write("Start thread #" + str(proc) + " " + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())) + '\n' + "-------------------------------------------------\n")
            for i in partition.index:
                try:
                    l.write("Convert record #" + str(i) + " " + str(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())) + "\n")
                    if isCorrectSequence(partition.loc[i, sequence_column], alphabets[alphabet]):

                        if alphabet == "protein":
                            flavor = 0
                        elif alphabet == "DNA":
                            flavor = 2
                        else:
                            flavor = 6

                        molecule = Chem.MolFromSequence(partition.loc[i, sequence_column],flavor = flavor)
                        molecule.SetProp('_Name', partition.loc[i, sequence_column])
                        molecule.SetProp('MolFileInfo',"Developed by S.I. Zhuravleva.")
                        #molecule.SetProp('Peptide', partition.loc[i, sequence_column])
                        if isCharged and alphabet == "protein":
                            chargePeptide(molecule)

                        # o.write(Chem.MolToMolBlock(molecule, forceV3000 = True)+"\n")
                        o.write(Chem.MolToMolBlock(molecule, forceV3000=True))
                        for c in partition.columns:
                            if c == sequence_column:
                                o.write(">  <" + "SEQUENCE" + ">\n")
                                o.write(str(partition.loc[i, c]) + "\n\n")
                                continue
                            if c == mutation_column:
                                o.write(">  <" + c + ">\n")
                                o.write(str(partition.loc[i, c]) + "\n\n")
                                o.write(">  <" + "POSITION" + ">\n")
                                pos = str(partition.loc[i, c])[3:-1]
                                o.write(pos + "\n\n")
                                continue
                            o.write(">  <" + c + ">\n")
                            o.write(str(partition.loc[i, c]) + "\n\n")
                        o.write("$$$$\n")
                        l.write("SUCCESS!\n")
                    else:
                        l.write("FAILED!\n")
                        f.write(str(i) + "\t" + partition.loc[i, sequence_column] + "\n")
                except Exception as er:
                    l.write("Exception on record #" + str(i) + "\t" + partition.loc[i, sequence_column] + "\n")
                    traceback.print_exc()
                    l.write(str(er) + "\n")
                    f.write(str(i) + "\t" + partition.loc[i, sequence_column] + "\n")
                finally:
                    l.write("-------------------------------------------------\n")
    except BaseException as e:
        print("Something went wrong in thread #" + str(proc)+"\n")
        traceback.print_exc()
        sys.exit(1)


def start_fun(configuration):
    start_time = time.time()

    if os.path.exists(configuration):
        try:
            with open(configuration, "r") as cfg:

                '''
                    Parsing arguments
                '''
                parameters = json.loads(''.join(cfg.readlines()))

                input = parameters["input"]
                output = parameters["output"]
                sequence_column = parameters["column_seq"]
                mutation_column = parameters["column_mut"]

                if "charged" in parameters.keys():
                    isCharged = parameters["charged"]
                else:
                    isCharged = False

                if "alphabet" in parameters.keys():
                    alphabet = parameters["alphabet"]
                else:
                    alphabet = "protein"

                if "threads" in parameters.keys():
                    number_threads = parameters["threads"]
                else:
                    number_threads = 1

                if "separator" in parameters.keys():
                    separator = parameters["separator"]
                else:
                    separator = ";"

                if "filename" in parameters.keys():
                    output_filename = parameters["filename"]
                else:
                    output_filename = os.path.splitext(os.path.basename(input))[0]

                '''
                      Validate arguments
                '''

                if not os.path.exists(input) or not os.path.isfile(input):
                    raise BaseException("Input file doesn`t exist or it isn`t file.")
                else:
                    table = pd.read_csv(input, sep=separator, header=0)
                    if sequence_column not in table.columns:
                        raise BaseException("Table doesn`t contain such column.")

                if not os.path.exists(output):
                    raise BaseException("Output directory doesn`t exist.")

                if type(isCharged) != bool:
                    raise BaseException("Which type amino acids residues it should use?")

                if alphabet not in alphabets.keys():
                    raise BaseException("Invalid alphabet. Allow only " + str(list(alphabets.keys())))

                if number_threads < 1 or number_threads > 2 * multiprocessing.cpu_count():
                    raise BaseException("Too many threads. Please use from 1 to " + str(2 * multiprocessing.cpu_count()) + " threads.")


            total = len(table.index)
            size_part = total // number_threads
            size_last_part = size_part + (total - size_part * number_threads)

            # procs - количество ядер
            # calc - количество операций на ядро

            processes = []

            # делим вычисления на количество ядер
            for proc, start in zip(range(number_threads), range(0, total, size_part)):
                if proc == number_threads - 1:
                    partition = table[start:start + size_last_part]
                else:
                    partition = table[start:start + size_part]

                p = multiprocessing.Process(target=write, args=(proc, partition, output, sequence_column, mutation_column, isCharged, alphabet, output_filename))
                processes.append(p)
                p.start()

            # Ждем, пока все ядра
            # завершат свою работу.
            for p in processes:
                p.join()

            # Merge all
            out_logs = output + "/logs/"

            sdf = os.path.join(output, output_filename + ".sdf")
            total_log = os.path.join(out_logs, output_filename + "_log" + ".txt")
            total_fail = os.path.join(out_logs, output_filename + "_failed" + ".txt")
            with open(sdf, "w", encoding="utf-8") as final, open(total_log, "w", encoding="utf-8") as log, open(total_fail, "w",encoding="utf-8") as failed:
                for i in range(number_threads):
                    out_t = os.path.join(out_logs, output_filename + "_thread_" + str(i) + ".sdf")
                    #os.mkdir(out_logs + output_filename + "_thread_" + str(i) + ".sdf")
                    log_t = os.path.join(out_logs, output_filename + "_thread_" + str(i) + "_log" + ".txt")
                    #os.mkdir(out_logs + output_filename + "_thread_" + str(i) + "_log" + ".txt")
                    fail_t = os.path.join(out_logs, output_filename + "_thread_" + str(i) + "_failed" + ".txt")
                    #os.mkdir(out_logs + output_filename + "_thread_" + str(i) + "_failed" + ".txt")
                    with open(out_t, "r", encoding="utf-8") as o, open(log_t, "r", encoding="utf-8") as l, open(fail_t, "r", encoding="utf-8") as f:
                        final.write("".join(o.readlines()))
                        log.write("".join(l.readlines()))
                        failed.write("".join(f.readlines()))
            print("Success")
        except BaseException as e:
            print("Something went wrong\n")
            traceback.print_exc()
            sys.exit(1)
        finally:
            end_time = time.time()
            print("--- %s seconds ---" % (end_time - start_time))
    else:
        print("Config doesn`t exist.")
        end_time = time.time()
        print("--- %s seconds ---" % (end_time - start_time))
        sys.exit(1)


if __name__ == '__main__':

    '''
        Path to config file
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Path to config file.")
    args = parser.parse_args()

    config = args.config
    start_fun(config)

