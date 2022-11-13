import pandas as pd
import os
import glob


def decisionTree(df, histone):
    nbColumns = ['Sequest q-value', 'Byonic Score', 'Byonic q-value']
    for column in df.columns.tolist():
        if column in nbColumns:
            df[column] = df[column].fillna(0)
        else:
            df[column] = df[column].fillna("")

    label = []
    for index, row in df.iterrows():
        if row["Proteins"] != 1:
            label.append("Rejected")
        elif row["Sequest Score"].upper() == "HIGH":
            label.append("Conserved")
        elif row["Byonic Score"] >= 300:
            label.append("Conserved")
        elif row["Sequest q-value"] <= 0.05 and row["Byonic q-value"] <= 0.05:
            label.append("Conserved")
        else:
            label.append("Verification") 
    df["Label"] = label
    if histone == "H4":
        df["Label"] = ["Conserved" for value in label]
    return df

def finalFilter(inDir):
    files = glob.glob(f"{inDir}/*analyzed.csv")
    for f in files:
        histone = os.path.basename(f).split("_")[1]
        print(histone)
        print(f)
        df = pd.read_csv(f, sep=";")
        df = decisionTree(df, histone)
        df.to_csv(f.replace("analyzed", "filtered"), sep=";", index=None)


if __name__ == "__main__":
    finalFilter("/home/kevin/Bureau/Aurelie/Donnees/PTM_figures/Mai_Grami/Delta")
    finalFilter("/home/kevin/Bureau/Aurelie/Donnees/PTM_figures/Mai_Grami/WT")
    finalFilter("/home/kevin/Bureau/Aurelie/Donnees/PTM_figures/Mai_Grami/OE")