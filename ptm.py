import pandas as pd
import re
import os
import glob
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
import subprocess
import copy

# matplotlib.rcParams.update({'font.size': 20})

def loadDF(path):
    return pd.read_excel(path)

def getMasterProtein(df, column, value):
    lastIndex = df.shape[0]
    df = df[df[column] == value]
    proteins = df["Accession"].tolist()
    mastersPos = df.index.values.tolist()
    indexTable = {"Accession": proteins, "Start":[], "End":[]}
    for i in range(len(mastersPos) - 1):
        indexTable["Start"].append(mastersPos[i] + 1)
        indexTable["End"].append(mastersPos[i+1])
    indexTable["Start"].append(mastersPos[-1])
    indexTable["End"].append(lastIndex)
    indexTable = pd.DataFrame.from_dict(indexTable)
    return indexTable

def countMasterProtein(inDir, organism):
    hcl = []
    h2so4 = []
    if organism == "Grami":
        conditions = getConditionsFile("/home/aurelie/Bureau/Aurelie/conditions.yaml")
    else:
        conditions = getConditionsFile("/home/aurelie/Bureau/Aurelie/fuji_conditions.yaml")
    for key in conditions.keys():
        if "HCl" in conditions[key]:
            hcl.append(key)
        elif "H2SO4" in conditions[key]:
            h2so4.append(key)
        else:
            pass
    hclDf = pd.DataFrame()
    h2so4Df = pd.DataFrame()
    files = sorted(glob.glob(f"{inDir}/*[1-9].xlsx"))
    for f in files:
        num = os.path.basename(f).split("_")[1].split(".")[0]
        df = loadDF(f)
        df = df.replace({"None":"Master Protein"})
        df = df.replace({"Master Protein Candidate":"Master Protein"})
        indexTable = getMasterProtein(df, "Master", 'Master Protein')
        if num in hcl:
            hclDf = pd.concat([hclDf, indexTable])
        elif num in h2so4:
            h2so4Df = pd.concat([h2so4Df, indexTable])
        else:
            pass
    print(hclDf)
    print(f"{organism} : HCl ({len(set(hclDf['Accession']))}) / H2SO4 ({len(set(h2so4Df['Accession']))})")

def parseMasterProtein(df, mastersPos):
    for i in range(len(mastersPos) - 1):
        accessionIndex = df.columns.tolist().index("Accession")
        start = mastersPos[i] + 1
        end = mastersPos[i+1]
        masterProt = df.iloc[start:end]
        return masterProt

def changeHeader(df):
    header = df.iloc[0]
    header = header.tolist()
    df2 = df.iloc[1:]
    df2.columns = header
    return df2

def getRelativePos(posListing):
    formatedPos = []
    i = 0
    while (i < len(posListing)):
        pos = ""
        while(posListing[i][-1] != "]"):
            pos += posListing[i]
            pos += ";"
            i += 1
        else:
            pos += posListing[i]
            formatedPos.append(pos)
        i += 1
    formatedPos = [re.search(r"\[.*\]", val).group() for val in formatedPos]
    formatedPos = ["".join(val[1:len(val)-1]).replace(" ", "") for val in formatedPos]
    return formatedPos


def getAbsolutePos(relativePosList, sequence, peptide):
    preAbsoluPos = sequence.find(peptide)
    absoluPosListing = []
    for pos in relativePosList:
        finalAbsoluPos = ""
        for p in pos.split(";"):
            if p == "N-Term":
                absoluPos = "N-Term"
            elif "/" in p:
                absoluPos = ""
                for val in p.split("/"):
                    absoluPos += val + "0"
                    absoluPos += "/"
                absoluPos = absoluPos[0:-1]
            else:
                aa = p[0]
                if len(p) == 1:
                    absoluPos = 0
                else:
                    absoluPos = preAbsoluPos + int(p[1:]) - 1
                finalAbsoluPos += aa
            finalAbsoluPos += str(absoluPos)
            finalAbsoluPos += ";"
        absoluPosListing.append(finalAbsoluPos[0:len(finalAbsoluPos)-1])
    return absoluPosListing

def replaceRelativeByAbosute(modifications, relativePosList, absoluPosList):
    modifications = modifications.split("; ")
    r = re.compile("[0-9]x*")
    modifications = list(filter(r.match, modifications))
    finalModif = ""
    for i in range(len(modifications)):
            finalModif += modifications[i].split(" ")[0]
            finalModif += " [{}]; ".format(absoluPosList[i])
    return finalModif[:-2]


def modificationListing(df, accession, dfGlobal, accessionData, condition, date, filename):
    sequence = accessionData["Sequence"]
    protein = accessionData["Protein"]
    dfGlobal = dfGlobal[dfGlobal["Accession"] != accession]
    data = {"Peptide": [], "Relative modification": [],"Absolute modification": [], "Shared count": []}
    for index, row in df.iterrows():
        peptide = row["Annotated Sequence"].split(".")[1]
        if row["Modifications"] == "":
            relativeModification = ""
            absoluteModification = ""
        else:
            relativeModification = row["Modifications"]
            ptms = relativeModification.split(";")
            relativePosList = getRelativePos(ptms)
            absolutePos = getAbsolutePos(relativePosList, sequence, peptide)
            absoluteModification = replaceRelativeByAbosute(relativeModification, relativePosList, absolutePos)
        data["Peptide"].append(peptide)
        data["Relative modification"].append(relativeModification)
        data["Absolute modification"].append(absoluteModification)
        count = dfGlobal[dfGlobal["Accession"] == row["Annotated Sequence"]].shape[0]
        data["Shared count"].append(count)
    final = pd.DataFrame.from_dict(data)
    final["Sequest Score"] = df["Confidence (by Search Engine): Sequest HT"].tolist()
    final["Proteins"] = df["# Proteins"].tolist()
    final["Sequest q-value"] = df["Percolator q-Value (by Search Engine): Sequest HT"]
    if "Byonic Score (by Search Engine): A5 PMI-Byonic" in df.columns.tolist():
        final["Byonic Score"] = df["Byonic Score (by Search Engine): A5 PMI-Byonic"].tolist()
        final["Byonic q-value"] = df["q-Value 2D (by Search Engine): A5 PMI-Byonic"]
    elif "Byonic Score (by Search Engine): A4 PMI-Byonic" in df.columns.tolist():
        final["Byonic Score"] = df["Byonic Score (by Search Engine): A4 PMI-Byonic"].tolist()
        final["Byonic q-value"] = df["q-Value 2D (by Search Engine): A4 PMI-Byonic"]
    else:
        final["Byonic Score"] = [""] * df.shape[0]
    final.insert(0, "Accession", [accession] * df.shape[0])
    final["Protein"] = [protein] * df.shape[0]
    final["Condition"] = [condition] * df.shape[0]
    final["Date"] = [date] * df.shape[0]
    final["Original file"] = [filename] * df.shape[0]
    return final


def parseFASTA(fastaFile):
    data = {}
    fasta = open(fastaFile, "r")
    lines = fasta.readlines()
    i = 0
    for i in range(len(lines)):
        if ">" in lines[i]:
            accession = lines[i].split("|")[1]
            protein = lines[i].split(" ")[2]
            data[accession] = {}
            data[accession]["Protein"] = protein
            sequence = ""
            j = i+1
            while (">" not in lines[j]) and (j < len(lines) - 1):
                sequence += lines[j].rstrip()
                j += 1
            data[accession]["Sequence"] = sequence
    sequence += lines[-1].rstrip()
    data[accession]["Sequence"] = sequence
    return data

def getConditionsFile(conditionsFile):
    conditionsFile = open(conditionsFile, "r")
    lines = conditionsFile.readlines()
    data = {}
    for l in lines:
        key = l.rstrip().split(" : ")[0]
        value = l.rstrip().split(" : ")[1]
        data[key] = value
    return data

def saveDifferentFiles(outDir, accessions, finalDF):
    for accession in accessions.keys():
        protein = accessions[accession]["Protein"]
        histoneFinalDF = finalDF[finalDF["Accession"] == accession]
        protein = protein.replace(".", "")
        filename = "{}/{}_{}_analyzed.csv".format(outDir, accession, protein)
        histoneFinalDF.to_csv(filename, sep=";", index=None)

def mergeConditions(accessions, outDir):
    for accession in accessions.keys():
        filename = "{}/{}_{}_analyzed.csv".format(outDir, accession, accessions[accession]["Protein"].replace(".", ""))
        histoneFinalDF = pd.read_csv(filename, delimiter=";")
        duplicate = histoneFinalDF.duplicated(subset=["Accession", "Peptide", "Relative modification", "Absolute modification"]).tolist()
        histoneFinalDF["Duplicate"] = duplicate
        newHistoneFinalDF = pd.DataFrame()
        sharedConditions = histoneFinalDF[histoneFinalDF["Duplicate"] == True]
        notSharedConditions = histoneFinalDF[histoneFinalDF["Duplicate"] == False]
        sharedConditions.fillna("", inplace=True)
        for index, row in sharedConditions.iterrows():
                duplicatedPeptide = sharedConditions[(sharedConditions["Accession"] == row["Accession"]) & (sharedConditions["Peptide"] == row["Peptide"]) & (sharedConditions["Relative modification"] == row["Relative modification"]) & (sharedConditions["Absolute modification"] == row["Absolute modification"])]
                conditions = set(duplicatedPeptide["Condition"].tolist())
                dates = set(duplicatedPeptide["Date"].tolist())
                conditions = ",".join(conditions)
                dates = ",".join(dates)
                duplicatedPeptide["Condition"] = [conditions] * duplicatedPeptide.shape[0]
                duplicatedPeptide["Date"] = [dates] * duplicatedPeptide.shape[0]
                newHistoneFinalDF = pd.concat([newHistoneFinalDF, duplicatedPeptide])
        newHistoneFinalDF = pd.concat([newHistoneFinalDF, notSharedConditions])
        newHistoneFinalDF.pop("Duplicate")
        newHistoneFinalDF.drop_duplicates(inplace=True)
        outFile = "{}/{}_{}_analyzed_merged.csv".format(outDir, accession, accessions[accession]["Protein"].replace(".", ""))
        newHistoneFinalDF.to_csv(outFile, sep=";", index=None)

def mapping(outDir, accessions, pattern):
    accessions = parseFASTA(accessions)
    text = ""
    if pattern == "filtered":
        filename = "_filtered.csv"
    else:
        filename = "_analyzed_merged.csv"
    for accession in accessions.keys():
        fileResults = "{}/{}_{}{}".format(outDir, accession, accessions[accession]["Protein"].replace(".", ""), filename)
        histoneFinalDF = pd.read_csv(fileResults, delimiter=";")
        if pattern == "filtered":
            histoneFinalDF = histoneFinalDF[histoneFinalDF["Label"].str.contains("Conserved")]
        if histoneFinalDF.shape[0] > 0:
            sequence = accessions[histoneFinalDF["Accession"].tolist()[0]]["Sequence"]
            text += sequence + "\n"
            sequenceIndex = [i for i in range(len(sequence))]
            sequenceMatch = []
            for index, row in histoneFinalDF.iterrows():
                match = sequence.find(row["Peptide"])
                sequenceMatch += [match + i for i in range(len(row["Peptide"]))]
            coverage = len(set(sequenceIndex).intersection(set(sequenceMatch))) / len(sequence) * 100
            message = "Couverture de {}% pour {}".format(coverage, accession)
            text += "################################" + "\n"
            text += message + "\n"
            text += "################################" + "\n"
        else:
            text += accessions[accession]["Sequence"] + "\n"
            coverage = 0
            message = "Couverture de {}% pour {}".format(coverage, accession)
            text += "################################" + "\n"
            text += message + "\n"
            text += "################################" + "\n"
    f = open(f"{outDir}/couverture.txt", "w")
    f.write(text)
    f.close()

def formateRawData(filename, accessions, condition, date):
    df = loadDF(filename)
    df = df.replace({"None":"Master Protein"})
    df = df.replace({"Master Protein Candidate":"Master Protein"})
    indexTable = getMasterProtein(df, "Master", 'Master Protein')
    intermediateDF = pd.DataFrame()
    for accession in accessions.keys():
        interestingAccession = indexTable[indexTable["Accession"] == accession]
        if interestingAccession.shape[0] != 0:
            start = interestingAccession["Start"].tolist()[0]
            end = interestingAccession["End"].tolist()[0]
            peptidesDF = df.iloc[start:end]
            peptidesDF = peptidesDF.dropna(axis=1, how="all")
            peptidesDF = changeHeader(peptidesDF)
            peptidesDF.fillna("", inplace=True)
            accessionData = accessions[accession]
            histoneDF = modificationListing(peptidesDF, accession, df, accessionData, condition, date, filename)
            intermediateDF = pd.concat([intermediateDF, histoneDF])
    return intermediateDF

def main(dataPath, fastaPath, conditionsPath, month):
    finalDates = {"Mars":"03.03.2020", "Mai":"07.05.2021", "Decembre":"16.12.2020", "Septembre":"30.09.2021"}
    outDir = os.path.dirname(dataPath)
    files = glob.glob(dataPath)
    accessions = parseFASTA(fastaPath)
    conditions = getConditionsFile(conditionsPath)
    finalDF = pd.DataFrame()
    for f in files:
        if len(f.split("_")) == 2:
            condition = f.split("_")[1]
            condition = condition.split(".x")[0]
            finalCondition = conditions[condition]
            finalDate = finalDates[month]
            df = formateRawData(f, accessions, finalCondition, finalDate)
        else:
            condition = os.path.basename(f).split(".")[0].split("_")[1]
            finalCondition = conditions[condition]
            finalDate = finalDates[month]
            df = formateRawData(f, accessions, finalCondition, finalDate)
            finalDF = pd.concat([finalDF, df])
    saveDifferentFiles(outDir, accessions, finalDF)
    mergeConditions(accessions, outDir)

def compareDF(df1Path, df2Path, comparisonType):
    df1 = pd.read_csv(df1Path, delimiter=";")
    df2 = pd.read_csv(df2Path, delimiter=";")
    df = df1.merge(df2, how="outer", on=["Accession", "Peptide", "Relative modification", "Absolute modification"], indicator=True)
    outFile = "{}_{}.csv".format("_".join(df1Path.split("_")[:2]), comparisonType)
    if comparisonType == "intersection":
        df = df[df["_merge"] == "both"]
    elif comparisonType == "left":
        df = df[df["_merge"] == "left_only"]
    elif comparisonType == "right":
        df = df[df["_merge"] == "right_only"]
    else:
        print("Please choose valid option (intersection, left, right)")
    df.to_csv(outFile, sep=";", index=None)

def mainCompare(dir1, dir2):
    histones = ["A0A098DX13_H4", "I1S219_H1", "Q4HTT1_H2A", "Q4HTT2_H2B", "Q4IER8_H3", "Q4IMD1_H2AZ", "V6R3L4_H4"]
    for histone in histones:
        f1 = glob.glob("{}{}*.csv".format(dir1, histone))[0]
        f2 = glob.glob("{}{}*.csv".format(dir2, histone))[0]
        compareDF(f1, f2, "intersection")
        compareDF(f1, f2, "left")
        compareDF(f1, f2, "right")

def splitPerCondition(df, typeCondition):
    if typeCondition == "extraction":
        return ["HCl", "H2SO4"]
    else:
        conditions = df["Condition"].str.split(",", expand=True)
        conditionsList = []
        for col in conditions.columns.tolist():
            conditionsList += conditions[col].tolist()
        conditionsList = [condition for condition in conditionsList if condition != None]
        conditionsList = list(set(conditionsList))
        if typeCondition == "splitted":
            return conditionsList
        else:
            conditionsList = [condition.split(" ")[1] for condition in conditionsList]
            if typeCondition == "enzyme":
                return list(set(conditionsList))
            else:
                conditionsList = list(set(conditionsList))
                concatCond = [condition for condition in conditionsList if condition != "Trypsin"]
                conditionsList = ["|".join(concatCond), "Trypsin"]
                return conditionsList

def mappingStats(fastaFile, resFile, conditionType, yMax, outDir):
    if "I1S219" in resFile:
        plt.subplots(figsize=(25, 3))
    else:
        plt.subplots(figsize=(20, 3))
    matplotlib.rc('font',family='Courier')
    fasta = parseFASTA(fastaFile)
    id = os.path.basename(resFile).split("_")[0]
    sequence = fasta[id]["Sequence"]
    coverage = []
    condData = []
    sequenceCol = []
    final = pd.DataFrame()
    if id in resFile:
        df = pd.read_csv(resFile, sep=";")
        df = df[df["Label"].str.contains("Conserved")]
        df["Condition"] = df["Condition"].str.replace("trypsin", "Trypsin")
        sequenceCount = {}
        if conditionType != "all":
            conditions = splitPerCondition(df, conditionType)
            if conditionType == "trypsin":
                boxplotPeptideLength(outDir, df, conditions)
            for condition in conditions:
                for i in range(len(sequence)):
                    sequenceCount[i] = 0
                cond = df[df["Condition"].str.contains(condition)]
                peptides = cond["Peptide"].tolist()
                for peptide in peptides:
                    match = sequence.index(peptide)
                    for aa in peptide:
                        sequenceCount[match] += 1
                        match += 1
                coverage += list(sequenceCount.values())
                condData += [condition] * len(list(sequenceCount.values()))
                sequenceCol += list(sequence)
                final[condition] = list(sequenceCount.values())
        else:
            for i in range(len(sequence)):
                sequenceCount[i] = 0
            peptides = df["Peptide"].tolist()
            for peptide in peptides:
                match = sequence.index(peptide)
                for aa in peptide:
                    sequenceCount[match] += 1
                    match += 1
            coverage += list(sequenceCount.values())
            condData += ["All"] * len(list(sequenceCount.values()))
            sequenceCol += list(sequence)
            final["All"] = list(sequenceCount.values())
    a = [final[col].tolist() for col in final.columns.tolist()]
    labels = sorted(final.columns.tolist())
    a.reverse()
    if len(a) > 0:
        if "I1S219" in resFile:
            size = 10
        else:
            size = 12
        h = sns.heatmap(a, xticklabels=list(sequence), linewidths=0.5, yticklabels=labels, cmap="Reds")
        h.tick_params(axis='x', which='major', labelsize=size, labelrotation=0)
        h.tick_params(axis='y', which='major', labelsize=size, labelrotation=0)
        h.figure.savefig(f'{outDir}/{id}_{conditionType}.png')
        plt.clf()
    else:
        print(f"#######################################NO COVERAGE FOR {resFile}")

def boxplotPeptideLength(outDir, df, conditions):
    final = pd.DataFrame()
    length = []
    conditionsRenamed = []
    for condition in conditions:
        cond = df[df["Condition"].str.contains(condition)]
        length += [len(sequence) for sequence in cond["Peptide"].tolist()]
        conditionsRenamed += [condition] * cond.shape[0]
    final["Condition"] = conditionsRenamed
    final["Length"] = length
    final.to_csv(f"{outDir}/length.tsv", sep="\t", index=None)
    subprocess.call(f"Rscript /home/aurelie/Bureau/Aurelie/boxplot.R {outDir}/length.tsv {outDir}", shell=True)

def sortDf(df, organism):
    gramiName = {"Q4HTT1" : "H2A", "Q4IMD1" : "H2A.Z", "Q4HTT2" : "H2B", "Q4IER8": "H3", "V6R3L4" : "H4.1", "A0A098DX13" : "H4.2", "I1S219" : "H1"}
    if "Fujikuroi" in organism:
        fuji = {"H2A" : 1, "H2A.Z" : 5, "H2B" : 2, "H3" : 3, "H4.1" : 4, "H1" : 0}
        df[1] = df[1].str.replace("H2AZ", "H2A.Z")
        grami = {"H2A" : 1, "H2A.Z" : 6, "H2B" : 2, "H3" : 3, "H4.1" : 4, "H4.2" : 5, "H1" : 0}
        df["sort"] = [fuji[value] for value in df[1].tolist()]
    else:
        grami = {"H2A" : 1, "H2A.Z" : 6, "H2B" : 2, "H3" : 3, "H4.1" : 4, "H4.2" : 5, "H1" : 0}
        if df[1].tolist()[0] not in list(gramiName.keys()):
            df[1] = df[1].str.replace("H2AZ", "H2A.Z")
            df["sort"] = [grami[value] for value in df[1].tolist()]
        else:
            df[1] = [gramiName[value] for value in df[1].tolist()]
            df["sort"] = [grami[value] for value in df[1].tolist()]
    df = df.sort_values(by=["sort"])
    df.pop("sort")
    return df

def barplotCoverage(coverageFile, outDir, organism):
    plt.subplots(figsize=(8.4, 6.8))
    subprocess.call(f"grep '[0-9]*\.*[0-9]*%' {coverageFile} | cut -d ' ' -f 3,5 > cov.txt", shell=True)
    df = pd.read_csv("cov.txt", sep=" ", header=None)
    df[0] = df[0].str.split("%", expand=True)[0].tolist()
    df[0] = df[0].astype(float)
    df[0] = [round(value, 1) for value in df[0].tolist()]
    df = sortDf(df, organism)
    df.columns = ["Coverage (%)", "Histone"]
    ax = sns.barplot(x="Histone", y="Coverage (%)", data=df, palette="Set3")
    ax.set_xlabel("Histone",fontsize=20)
    ax.set_ylabel("Coverage (%)",fontsize=20)
    ax.tick_params(labelsize=18)
    plt.savefig(f'{outDir}/coverage.png')
    plt.clf()

def cleanModifications(finalModifications, histone, inDir):
    dico = defaultdict(list)
    finalModifications = set(finalModifications)
    finalModifications = [value for value in finalModifications if not "Met-loss" in value]
    newFinalModifications = []
    for value in finalModifications:
        if re.findall('[A-Z]0', value):
            if "Oxidation" in value:
                nb, o, modif = value.split("x")
                nb = int(nb)
                nb -= 1
                newFinalModifications.append(f"{nb}x{o}x{modif}")
            else:
                nb, modif = value.split("x")
                nb = int(nb)
                nb -= 1
                newFinalModifications.append(f"{nb}x{modif}")
        else:
            newFinalModifications.append(value)
    for modif in newFinalModifications:
        if "Oxidation" in modif:
            nb, o, modif = modif.split("x")
            modif, aa = modif.split(" ")
            modif = f"{o}x{modif}"
            aa = aa[1:-1]
            aa = aa.split(";")
        else:
            nb, modif = modif.split("x")
            modif, aa = modif.split(" ")
            aa = aa[1:-1]
            aa = aa.split(";")
        for a in aa:
            dico[modif].append(a)
    dico2 = copy.deepcopy(dico)
    for key in dico2.keys():
        dico2[key] = [value for value in dico2[key] if not re.findall('[A-Z]0', value)]
        dico2[key] = [sorted(list(set(dico2[key])))]
    dicoToDf(dico2, histone, inDir)
    for key in dico.keys():
        dico[key] = [value for value in dico[key] if not re.findall('[A-Z]0', value)]
        dico[key] = [value for value in dico[key] if not "N-Term" in value]
        dico[key] = [len(list(set(dico[key])))]
    df = pd.DataFrame().from_dict(dico, orient = "index")
    if df.shape[0] != 0:
        df.columns = ["Count"]
        df["Modif"] = df.index
        df = df[["Modif", "Count"]]
    return df

def pieChartModif(resFile, inDir):
    plt.subplots(figsize=(6.4, 4.8))
    colorsDict = colorPalette()
    id = os.path.basename(resFile).split("_")[0]
    df = pd.read_csv(resFile, sep=";")
    df = df[df["Label"].str.contains("Conserved")]
    if df.shape[0] > 0 and df[~pd.isna(df["Absolute modification"])].shape[0] > 0:
        df = df[~pd.isna(df["Absolute modification"])]
        modifications = list(set(df["Absolute modification"]))
        finalModifications = []
        for modification in modifications:
            if len(modification.split("; ")) > 1:
                for value in modification.split("; "):
                    finalModifications.append(value)
            else:
                finalModifications.append(modification)
        finalModifications = list(set(finalModifications))
        histone = os.path.basename(resFile).split("_")[1]
        if "A0A" in resFile:
            histone = "H4.2"
        if "V6R3L4" in resFile:
            histone = "H4.1"
        df = cleanModifications(finalModifications, histone, inDir)
        if df.shape[0] != 0:
            df["Modif"] = df["Modif"].str.replace("GGQ", "Ubiquitylation")
            df["Modif"] = df["Modif"].str.replace("GlyGlyGln", "Ubiquitylation")
            df["Modif"] = df["Modif"].str.replace("Acetyl", "Acetylation")
            df["Modif"] = df["Modif"].str.replace("Formyl", "Formylation")
            df["Modif"] = df["Modif"].str.replace("Methyl", "Monomethylation")
            df["Modif"] = df["Modif"].str.replace("Dimethyl", "Dimethylation")
            df["Modif"] = df["Modif"].str.replace("Trimethyl", "Trimethylation")
            df["Modif"] = df["Modif"].str.replace("Crotonyl", "Crotonylation")
            df["Modif"] = df["Modif"].str.replace("Phospho", "Phosphorylation")
            df["Histone"] = [histone] * df.shape[0]
            df = df[df["Modif"] != "Oxidation"]
            pie = plt.pie(df["Count"].tolist(), colors = [colorsDict[modif] for modif in df["Modif"].tolist()], autopct='%.0f%%', textprops={'fontsize': 20})
            plt.title(os.path.basename(resFile).split("_")[1])
            plt.savefig(f'{os.path.dirname(resFile)}/{id}_pie.png')
            plt.clf()
        return df
    else:
        print(f"#######################################NO COVERAGE FOR {resFile}")
        return pd.DataFrame()

def pieChartHistone(df, outDir, organism):
    plt.subplots(figsize=(6.4, 4.8))
    count = []
    histones = set(df["Histone"])
    label = []
    for histone in histones:
        subDf = df[df["Histone"] == histone]
        label.append(histone)
        count.append(subDf["Count"].sum())
    histones = df[["Histone"]].drop_duplicates()
    histones.columns = [1]
    histones = sortDf(histones, organism)
    colors = sns.color_palette('Set3')[0:len(df["Modif"].tolist())]
    print("LABEL : ", label)
    print("COUNT : ", count)
    plt.pie(count, labels = label, colors = colors, autopct='%.0f%%', textprops={'fontsize': 20})
    plt.savefig(f'{outDir}/Histone_pie.png')
    plt.clf()

def computeStats(inDir, histones, legend, size, typeAnalysis):
    conditions = ["all", "extraction", "splitted", "enzyme", "trypsin"]
    gramiFiles = glob.glob(f"{inDir}/*filtered.csv")
    df = pd.DataFrame()
    for gramiFile in gramiFiles:
        for condition in conditions:
            mappingStats(histones, gramiFile, condition, 10, os.path.dirname(gramiFile))
    df = pd.DataFrame()
    for gramiFile in gramiFiles:
        df = pd.concat([df, pieChartModif(gramiFile, inDir)])
        df.to_csv(f"{inDir}/nested.tsv", sep="\t", index=None)
    mapping(inDir, histones, "filtered")
    pieChartHistone(df, inDir, histones)
    barplotCoverage(f"{inDir}/couverture.txt", inDir, histones)
    colocalizate(inDir)
    mergeColocalizate(inDir)
    heatmapColocalizate(inDir, histones)

def nestedPieBothOrganism(gramiDir, fujiDir):
    nestedGrami = pd.read_csv(glob.glob(f"{gramiDir}/nested.tsv")[0], sep="\t")
    nestedFuji = pd.read_csv(glob.glob(f"{fujiDir}/nested.tsv")[0], sep="\t")
    nestedGrami["Species"] = ["Fusarium graminearum"] * nestedGrami.shape[0]
    nestedFuji["Species"] = ["Fusarium fujikuroi"] * nestedFuji.shape[0]
    nested = pd.concat([nestedGrami, nestedFuji])
    nested.to_csv(f"{gramiDir}/nestedConcat.tsv", sep="\t", index=None)

def flatList(listValues):
    flat_list = []
    for listValue in listValues:
        for value in listValue:
            flat_list.append(value)
    return flat_list

def reshapeListofList(refList, flattedList):
    indexes = [0]
    reshapedList = []
    for liste in refList:
        indexes.append(indexes[-1] + len(liste))
    for i in range(len(indexes) - 1):
        reshapedList.append(flattedList[indexes[i]:indexes[i+1]])
    return reshapedList

def normValuesForNested(listValues, sumValues):
    for listValue in listValues:
        for i in range(len(listValue)):
            listValue[i] = listValue[i] / sum(sumValues) * 2 * np.pi
    return listValues

def getBarEdgesOrdinates(flattedList):
    flattedList = flattedList[:-1]
    ordinates = [0, flattedList[0]]
    for i in range(1, len(flattedList)):
        ordinates.append(flattedList[i] + ordinates[-1])
    return ordinates

def colorPalette():
    cmap1 = matplotlib.cm.get_cmap('Pastel1')
    cmap2 = matplotlib.cm.get_cmap('Pastel2')
    cmap3 = matplotlib.cm.get_cmap('tab10')
    cmap4 = matplotlib.cm.get_cmap('Accent')
    pastel1 = [matplotlib.colors.rgb2hex(cmap1(i)[:3]) for i in range(cmap1.N)]
    pastel2 = [matplotlib.colors.rgb2hex(cmap2(i)[:3]) for i in range(cmap2.N)]
    tab10 = [matplotlib.colors.rgb2hex(cmap3(i)[:3]) for i in range(cmap3.N)]
    accent = [matplotlib.colors.rgb2hex(cmap4(i)[:3]) for i in range(cmap4.N)]
    # modifs = ["Acetylation", "Crotonylation" ,"Dimethylation", "Trimethylation", "Methylation", "Phosphorylation", "Oxidation", "Formylation", "Ubiquitylation", "Met-loss", "Met-loss+Acetylation"]
    modifs = ["Acetylation", "Crotonylation" ,"Dimethylation", "Trimethylation", "Monomethylation", "Phosphorylation", "Formylation", "Ubiquitylation", "Met-loss", "Met-loss+Acetylation"]
    # colors = [pastel1[1], pastel1[6], pastel1[2], pastel2[4], pastel2[0], pastel1[3], pastel1[8], pastel1[7], pastel1[4], pastel1[0], pastel2[2]]
    colors = [pastel1[1], pastel1[6], pastel1[2], pastel2[4], pastel2[0], pastel1[3], pastel1[7], pastel1[4], pastel1[0], pastel2[2]]
    colors = [tab10[0], tab10[5], tab10[8], accent[0], tab10[2], tab10[4], tab10[6], tab10[1], pastel1[0], pastel2[2]] 
    return dict(zip(modifs, colors))
    #bleu : acetylation
    #marron : crotonylation
    #vert pastel1 : dimethylation
    #vert pastel2 1 : trimethylation
    #vert pastel2 2 : methylation
    #violet : phosporylation
    #gris : oxydation
    #rose : formylation
    #orange : ubiquitylation
    #rouge : met-loss

def colocalizate(inDir):
    tablesPath = glob.glob(f"{inDir}/*filtered.csv")
    for tablePath in tablesPath:
        df = pd.read_csv(tablePath, sep=";")
        df = df[df["Label"].str.contains("Conserved")]
        df = df[~pd.isna(df["Absolute modification"])]
        if df.shape[0] > 0:
            coloc1 = df[df["Absolute modification"].str.contains("; ")]
            coloc2 = df[~df["Absolute modification"].str.contains("1x")]
            colocalizated = pd.concat([coloc1, coloc2])
            tablePath = tablePath.replace("filtered", "coexisting")
            colocalizated.to_csv(tablePath, sep=";", index=None)

def mergeColocalizate(inDir):
    final = pd.DataFrame()
    subprocess.call(f"rm {inDir}/coexisting.csv", shell=True)
    tablesPath = glob.glob(f"{inDir}/*coexisting*")
    for tablePath in tablesPath:
        df = pd.read_csv(tablePath, sep=";")
        final = pd.concat([final, df])
    final.to_csv(f"{inDir}/coexisting.csv", sep=";", index=None)

def heatmapColocalizate(inDir, fasta):
    if not "Fujikuroi" in fasta:
        organism = "Grami"
    else:
        organism = "Fuji"
    dico = {"Grami": {"Q4HTT2" : (14, 12), "A0A098DX13" : (14, 7), "I1S219" : (20, 3), "Q4HTT1" : (12, 6), "Q4IER8" : (15, 14), "Q4IMD1" : (20, 3), "V6R3L4" : (10,5)}, \
        "Fuji" : {"H2B" : (12, 6), "H4" : (14, 3), "I1S219" : (20, 3), "H2A" : (13, 2), "H3" : (12, 4), "H2AZ" : (20, 3), "H1" : (20, 2)}}
    matplotlib.rc('font',family='Courier')
    coexisting = pd.read_csv(glob.glob(f"{inDir}/coexisting.csv")[0], sep=";")
    accessions = parseFASTA(fasta)
    for histone in set(coexisting["Accession"]):
        plt.subplots(figsize = dico[organism][histone])
        peptides = []
        subDf = coexisting[coexisting["Accession"] == histone]
        subDf = subDf.sort_values(by=["Peptide"])
        for peptide in subDf["Peptide"].tolist():
            peptides.append(mapColocalizate(peptide, accessions[histone]["Sequence"]))
        scores = setScoreFromPeptide(peptides)
        scores = setScoreForModif(subDf, scores)
        indexToRemove = filterOnlyOneModif(scores)
        colors = setColorForHeatmapColocalizate(scores)
        scores = attributeColor(scores, colors)
        indexToRemove = filterOnlyOneModif(scores)
        a = pd.DataFrame(scores)
        a = a[~a.index.isin(indexToRemove)]
        b = pd.DataFrame(peptides)
        b = b[~b.index.isin(indexToRemove)]
        a = a.drop_duplicates()
        b = b[b.index.isin(a.index.tolist())]
        dico[histone] = countDuplicatedPeptides(a)
        if a.shape[0] > 0:
            ax = sns.heatmap(a, annot=b, fmt="", cmap=colors, xticklabels=list(accessions[histone]["Sequence"]), yticklabels=[], vmin=0, cbar = False, annot_kws={"size": 10})
            ax.tick_params(axis='x', which='major', labelsize=12, labelrotation=0)
            plt.title(f"{histone}")
            plt.savefig(f'{inDir}/{histone}_colocalizated.png')
            plt.clf()

def countDuplicatedPeptides(df):
    dicoCount = defaultdict(lambda : 0)
    for i in range(df.shape[0]):
        scores = [str(val) for val in df.iloc[i].tolist()]
        scores = str(scores)
        dicoCount[scores] += 1
    return dicoCount

def attributeColor(scores, colors):
    indexColors = [i for i in range(len(colors))]
    values = []
    for score in scores:
        for value in score:
            values.append(value)
    values = list(set(sorted(values)))
    correspondance = dict(zip(values, indexColors))
    for score in scores:
        for i in range(len(score)):
            score[i] = correspondance[score[i]]
    return scores

def filterOnlyOneModif(scores):
    indexToRemove = []
    for i in range(len(scores)):
        if len(set(scores[i])) == 2:
            indexToRemove.append(i)
        elif len(set(scores[i])) == 3:
            value = [x for x in list(set(scores[i])) if x not in [0,1]]
            cptModif = 0
            for j in range(len(scores[i])):
                if scores[i][j] == value[0]:
                    cptModif += 1
            if cptModif <= 1:
                indexToRemove.append(i)
        else:
            pass
    return indexToRemove

def setColorForHeatmapColocalizate(scores):
    modifsColor = colorPalette()
    modifScore = {"Acetylation" : 2, "Crotonylation" : 3, "Dimethylation": 4, "Trimethylation" : 5, "Monomethylation" : 6, \
        "Phosphorylation": 7, "Formylation" : 9, "Ubiquitylation" : 10}
    modifScore = dict(zip(list(modifScore.values()), list(modifScore.keys())))
    colors = []
    values = []
    for score in scores:
        for value in score:
            values.append(value)
    values = sorted(list(set(values)))
    for value in values:
        if value == 0:
            colors.append("#ffffff")
        elif value == 1:
            colors.append("#ffffff")
        else:
            colors.append(modifsColor[modifScore[value]])
    return colors


def mapColocalizate(peptide, sequenceProtein):
    start = sequenceProtein.index(peptide)
    peptide = list(peptide)
    if start == 0:
        peptide += [" "] * (len(sequenceProtein) - len(peptide))
    else:
        if start + len(peptide) == len(sequenceProtein):
            peptide = [" "] * (len(sequenceProtein[:start]) + peptide)
        else:
            startMapping = [" "] * len(sequenceProtein[:start])
            endMapping = [" " ] * (len(sequenceProtein) - (len(peptide) + start))
            peptide = startMapping + peptide + endMapping
    return peptide

def setScoreFromPeptide(peptides):
    score = {" " : 0, 'C': 1, 'F': 1, 'T': 1, 'M': 1, 'E': 1, 'S': 1, 'A': 1, 'I': 1, \
        'R': 1, 'V': 1, 'N': 1, 'P': 1, 'W': 1, 'G': 1, 'K': 1, 'Y': 1, 'L': 1, 'Q': 1, 'D': 1, 'H': 1}
    scores = []
    for peptide in peptides:
        pepScore = []
        for aa in peptide:
            pepScore.append(score[aa])
        scores.append(pepScore)
    return scores

def setScoreForModif(df, scores):
    regexp = re.compile(r'[1-9][0-9]*')
    df = df.reset_index()
    modifScore = {"Acetylation" : 2, "Crotonylation" : 3, "Dimethylation": 4, "Trimethylation" : 5, "Monomethylation" : 6, \
        "Phosphorylation": 7, "Formylation" : 9, "Ubiquitylation" : 10}
    df["Absolute modification"] = df["Absolute modification"].str.replace("GGQ", "Ubiquitylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("GlyGlyGln", "Ubiquitylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Acetyl", "Acetylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Formyl", "Formylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Methyl", "Monomethylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Dimethyl", "Dimethylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Trimethyl", "Trimethylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Crotonyl", "Crotonylation")
    df["Absolute modification"] = df["Absolute modification"].str.replace("Phospho", "Phosphorylation")
    for index, row in df.iterrows():
        ptms = row["Absolute modification"].split(";")
        modifs = row["Absolute modification"].split("; ")
        ptmsPos = getRelativePos(ptms)
        for i in range(len(ptmsPos)):
            ptmPosSplitted = ptmsPos[i].split(";")
            for ptmPosSplit in ptmPosSplitted:
                if regexp.search(ptmPosSplit):
                    indexPos = int(regexp.search(ptmPosSplit).group(0))
                    for key in modifScore.keys():
                        if key in modifs[i]:
                            modif = key
                    scores[index][indexPos] = modifScore[modif]
                else:
                    pass
    return scores

def reshapeConditions(inDir):
    for f in glob.glob(f"{inDir}/*analyzed_merged.csv"):
        df = pd.read_csv(f, sep=";")
        df["Condition2"] = df["Condition"].str.split(",")
        df = df.explode("Condition2")
        df.pop("Condition")
        df.columns = [col if col != "Condition2" else "Condition" for col in df.columns]
        df.to_csv(f, sep=";", index=None)

def searchPeptideShared(inDir):
    files = glob.glob(f"{inDir}/*filtered.csv")
    for f in files:
        others = []
        df = pd.read_csv(f, sep=";")
        for index, row in df.iterrows():
            value = []
            if row["Proteins"] > 1:
                filename = row["Original file"]
                ref = pd.read_excel(filename)
                ref = ref.replace({"None":"Master Protein"})
                ref = ref.replace({"Master Protein Candidate":"Master Protein"})
                indexTable = getMasterProtein(ref, "Master", 'Master Protein')
                indexTable.to_csv("ok.tsv", sep="\t", index=None)
                proteins = ref["Master"].str.split(".", expand=True)[1].fillna("").tolist()
                indexToKeep = [i for i in range(len(proteins)) if proteins[i] == row["Peptide"]]
                for indexe in indexToKeep:
                    subIndexTable = indexTable[(indexe > indexTable["Start"]) & (indexe < indexTable["End"])]
                    value.append(subIndexTable["Accession"].tolist()[0])
                others.append("/".join(value))
            else:
                others.append("")
        df["Share_Protein"] = others
        df.to_csv(f, sep=";", index=None)

def concatMonths(fasta, inDir, species):
    folders = glob.glob(f"{inDir}/*{species}")
    dico = parseFASTA(fasta)
    accessions = dico.keys()
    subprocess.call(f"mkdir -p {inDir}/Final_Résultats/{species}/", shell=True)
    for accession in accessions:
        concat = pd.DataFrame()
        for folder in folders:
            f = glob.glob(f"{folder}/*{accession}*filtered.csv")
            if len(f) != 0:
                concat = pd.concat([concat, pd.read_csv(f[0], sep=";")])
        name = os.path.basename(f[0])
        concat.to_csv(f"{inDir}/{name}", sep=";", index=None)

def removeOxidation(inDir):
    files = glob.glob(f"{inDir}/*filtered.csv")
    for f in files:
        print(f)
        df = pd.read_csv(f, sep=";")
        for index, row in df.iterrows():
            for col in ["Absolute modification", "Relative modification"]:
                if not pd.isna(row[col]):
                    if "Oxidation" in row[col]:
                        modifs = row[col].split("; ")
                        if len(modifs) == 1:
                            df.at[index, col] = ""
                        else:
                            indexToPop = [i for i in range(len(modifs)) if "Oxidation" in modifs[i]]
                            modifs.pop(indexToPop[0])
                            df.at[index, col] = "; ".join(modifs)
        df.to_csv(f, sep=";", index=None)

def dicoToDf(dico, histone, inDir):
    dicoCleaned = {}
    for key in dico.keys():
        dicoCleaned[key] = dico[key][0]
    modif = []
    residue = []
    for key in dicoCleaned.keys():
        residue += dicoCleaned[key]
        modif += [key] * len(dicoCleaned[key])
    final = pd.DataFrame()
    final["Modification"] = modif
    final["Residue"] = residue
    final = final[final["Residue"] != "N-Term"]
    residues = final["Residue"].tolist()
    residues = [int(value[1:]) for value in residues]
    final["sort"] = residues
    final = final.sort_values(by = ["sort"])
    final = final[final["Modification"] != "Oxidation"]
    final.pop("sort")
    final["Modification"] = final["Modification"].str.replace("GGQ", "Ubiquitylation")
    final["Modification"] = final["Modification"].str.replace("GlyGlyGln", "Ubiquitylation")
    final["Modification"] = final["Modification"].str.replace("Acetyl", "Acetylation")
    final["Modification"] = final["Modification"].str.replace("Formyl", "Formylation")
    final["Modification"] = final["Modification"].str.replace("Methyl", "Monomethylation")
    final["Modification"] = final["Modification"].str.replace("Dimethyl", "Dimethylation")
    final["Modification"] = final["Modification"].str.replace("Trimethyl", "Trimethylation")
    final["Modification"] = final["Modification"].str.replace("Crotonyl", "Crotonylation")
    final["Modification"] = final["Modification"].str.replace("Phospho", "Phosphorylation")
    final["Histone"] = [histone] * final.shape[0]
    final["Reference"] = [""] * final.shape[0]
    final = final[["Histone", "Residue", "Modification", "Reference"]]
    final.to_csv(f"{inDir}/{histone}_DECOMPTE.tsv", sep = "\t", index = None)

def concatDf(inDir):
    files = glob.glob(f"{inDir}/*DECOMPTE.tsv")
    final = pd.DataFrame()
    for f in files:
        if not "H4_DECOMPTE" in f:
            final = pd.concat([final, pd.read_csv(f, sep = "\t")])
    final["Histone"] = final["Histone"].str.replace("H2AZ", "H2A.Z")
    sorter = ["H1", "H2A", "H2B", "H3", "H4.1", "H4.2", "H2A.Z"]
    final["Histone"] = final["Histone"].astype("category")
    final["Histone"] = final["Histone"].cat.set_categories(sorter)
    residues = final["Residue"].tolist()
    residues = [int(value[1:]) for value in residues]
    final["sort"] = residues
    final = final.sort_values(by=["Histone", "sort"])
    final.pop("sort")
    final = final.drop_duplicates(subset = final.columns.tolist())
    final.to_csv(f"{inDir}/FINAL_DECOMPTE.tsv", index = None, sep = "\t")

def filterCoExisting(inDir):
    files = glob.glob(f"{inDir}/*coexisting.csv")
    finalDf = pd.DataFrame()
    for f in files:
        tmpDf = pd.DataFrame()
        peptide = []
        histoneFinal = []
        modif = []
        if os.path.basename(f) != "coexisting.csv":
            histone = os.path.basename(f).split("_")[1]
            df = pd.read_csv(f, sep = "\t")
            if "A0A" in f:
                histone = "H4.2"
            if "V6R3L4" in f:
                histone = "H4.1"
            for index, row in df.iterrows():
                modifs = row["Absolute modification"].split("; ")
                finalModif = ""
                for value in modifs:
                    typeModif, residues = value.split(" ")
                    count = int(typeModif.split("x")[0])
                    noFinalResidues = len(residues.split(";"))
                    residues = residues[1:-1]
                    finalResidues = re.findall("[A-Z][1-9][0-9]*[0-9]*", residues)
                    finalModif += f"{len(finalResidues)}x{typeModif.split('x')[-1]} [{';'.join(finalResidues)}]; "
                keepModif = []
                for value in finalModif.split("; "):
                    toKeep = value.split("x")[-1].split(" ")[0]
                    if "idation" == toKeep:
                        pass
                    elif "Met-loss+Acetyl" == toKeep:
                        pass
                    elif "Met-loss" == toKeep:
                        pass
                    elif value == "":
                        pass
                    else:
                        keepModif.append(value)
                finalModif = "; ".join(keepModif)
                finalModif = [value for value in finalModif.split("; ") if not "O" in value]
                finalModif = "; ".join(finalModif)
                if finalModif != "":
                    total = sum([int(value.split("x")[0]) for value in finalModif.split("; ")])
                    if total > 1 or len(finalResidues) > 1:
                        peptide.append(row["Peptide"])
                        modif.append(finalModif)
            histoneFinal += [histone] * len(modif)
            tmpDf["Histone"] = histoneFinal
            tmpDf["Residue"] = modif
            tmpDf["Peptide"] = peptide
            finalDf = pd.concat([finalDf, tmpDf])
    finalDf["Residue"] = finalDf["Residue"].str.replace("GGQ", "ub")
    finalDf["Residue"] = finalDf["Residue"].str.replace("GlyGlyGln", "ub")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Acetyl", "ac")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Formyl", "formyl")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Methyl", "me1")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Dimethyl", "me2")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Trimethyl", "me3")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Crotonyl", "croto")
    finalDf["Residue"] = finalDf["Residue"].str.replace("Phospho", "ph")
    newResidue = []
    for residue in finalDf["Residue"].tolist():
        formated = []
        for value in residue.split("; "):
            modif, res = value.split(" ")
            modif = modif.split("x")[-1]
            res = res[1:-1]
            res = res.split(";")
            for val in res:
                if val != "":
                    formated.append(f"{val}{modif}")
        newResidue.append(" + ".join(formated))
    finalDf["New residue"] = newResidue
    finalDico = {}
    for histone in set(finalDf["Histone"]):
        subDf = finalDf[finalDf["Histone"] == histone]
        dico = defaultdict(int)
        for coloc in set(subDf["New residue"]):
            for value in subDf["New residue"].tolist():
                if coloc in value:
                    dico[coloc] += 1
        finalDico[histone] = dico
    dfCount = pd.DataFrame()
    for key in finalDico.keys():
        df = pd.DataFrame()
        df["Histone"] = [key] * len(list(finalDico[key].keys()))
        df["New residue"] = list(finalDico[key].keys())
        df["count"] = list(finalDico[key].values())
        dfCount = pd.concat([df, dfCount])
    finalDf = finalDf.merge(dfCount, on = ["Histone", "New residue"])
    finalDf = finalDf[["Histone", "New residue", "count"]]
    finalDf = finalDf.drop_duplicates()
    finalDf["Histone"] = finalDf["Histone"].str.replace("H2AZ", "H2A.Z")
    sorter = ["H1", "H2A", "H2B", "H3", "H4.1", "H4.2", "H2A.Z"]
    finalDf["Histone"] = finalDf["Histone"].astype("category")
    finalDf["Histone"] = finalDf["Histone"].cat.set_categories(sorter)
    finalDf = finalDf.sort_values(by = ["Histone"])
    finalDf.to_excel(f"{inDir}/coexisting.xlsx", index = None)#, sep = "\t", index = None)

def peptideLength(inDir, yaml):
    files = glob.glob(f"{inDir}/*filtered.csv")
    dico = {}
    lines = open(yaml, "r").readlines()
    final = pd.DataFrame()
    for line in lines:
        value = line.rstrip()
        dico[value.split(":")[0][:-1]] = value.split(":")[1][1:]
    for f in files:
        df = pd.read_csv(f, sep = ";")
        df["Original file"] = [os.path.basename(value).split("_")[1].split(".")[0] for value in df["Original file"].tolist()]
        df["Original file"] = [dico[value].replace("trypsin", "Trypsin") for value in df["Original file"].tolist()]
        final = pd.concat([final, df])
    df["Condition"] = df["Original file"].str.split(" ", expand = True)[1]
    df["Length"] = [len(value) for value in df["Peptide"].tolist()]
    df = df[["Condition", "Length"]]
    df.to_csv(f"{inDir}/peptideLength.tsv", sep = "\t", index = None)

if __name__ == "__main__":
    # main("./Donnees/PTM_figures/Mai_Fuji/*.xlsx", "Fujikuroi.fasta", "fuji_conditions.yaml", "Mai")
    # main("./Donnees/PTM_figures/Mai_Grami/*.xlsx", "Histones2.fasta", "conditions.yaml", "Mai")
    # main("./Donnees/PTM_figures/Mars_Grami/*.xlsx", "Histones2.fasta", "conditions.yaml", "Mars")
    # main("./Donnees/PTM_figures/Septembre_Grami/*.xlsx", "Histones2.fasta", "conditions.yaml", "Septembre")
    # main("./Donnees/PTM_figures/Decembre_Grami/*.xlsx", "Histones2.fasta", "conditions.yaml", "Decembre")
    # main("./Donnees/PTM_figures/Mai_Grami/Delta/*.xlsx", "Histones2.fasta", "conditions.yaml", "Mai")
    # main("./Donnees/PTM_figures/Mai_Grami/OE/*.xlsx", "Histones2.fasta", "conditions.yaml", "Mai")
    # main("./Donnees/PTM_figures/Mai_Grami/WT/*.xlsx", "Histones2.fasta", "conditions.yaml", "Mai")
    # main("./Donnees/Mai/*.xlsx", "Fujikuroi.fasta", "fuji_conditions.yaml")
    #main("./Donnees/Decembre/*.xlsx", "Histones.fasta", "conditions.yaml")

    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Mars_Grami", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Septembre_Grami", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Decembre_Grami", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Mai_Grami", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Decembre_Grami/OE", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Decembre_Grami/Delta", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Decembre_Grami/WT", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    computeStats("/home/aurelie/Bureau/Aurelie/Donnees/Concat_OE", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    computeStats("/home/aurelie/Bureau/Aurelie/Donnees/Concat_Delta", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    computeStats("/home/aurelie/Bureau/Aurelie/Donnees/Concat_WT", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/PTM_figures/Mai_Fuji", "/home/aurelie/Bureau/Aurelie/Fujikuroi.fasta", False, 0.3, 'Histone')
    # concatMonths("/home/aurelie/Bureau/Aurelie/Histones2.fasta", "/home/aurelie/Bureau/Aurelie/Donnees/Concat/", "Grami")
    # concatMonths("/home/aurelie/Bureau/Aurelie/Histones2.fasta", "/home/aurelie/Bureau/Aurelie/Donnees/Concat_OE", "Grami")
    # concatMonths("/home/aurelie/Bureau/Aurelie/Histones2.fasta", "/home/aurelie/Bureau/Aurelie/Donnees/Concat_WT", "Grami")
    # concatMonths("/home/aurelie/Bureau/Aurelie/Histones2.fasta", "/home/aurelie/Bureau/Aurelie/Donnees/Concat_Delta", "Grami")
    # concatMonths("/home/aurelie/Bureau/Aurelie/Histones2.fasta", "/home/aurelie/Bureau/Aurelie/Donnees/Concat_WithoutSeptembre/", "Grami")
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/Concat/", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')
    # computeStats("/home/aurelie/Bureau/Aurelie/Donnees/Concat_WithoutSeptembre/", "/home/aurelie/Bureau/Aurelie/Histones2.fasta", False, 0.3, 'Histone')