from src.entities import GraphStats, ImmuneNetwork
from src.creation.algorithms.simpleDistance import simpleDistance
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 
from src.analysis.visualization.multiGraphChart import multiGraphChart

class LatexPrettyPrinter:
    def generateLatexHeuristicHeader(self, results, heuristicStep=2, polish=True):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
        paramNamesDict = {
            "algorithm_name":"Nazwa algorytmu",
            "distance":"Funckja dystansu",
            "substitution_matrix":"Macierz podstawień",
            "threshold":"Próg podobieństwa"
        }
        
        if polish:
            paramNames = [ paramNamesDict[param] for param in paramNames]
        else:
            paramNames = [ name.replace("_", " ") for name in paramNames]
    
        headerPrefixEng = [
            "test group",
        ]
        headerPrefixPl = [
            "Grupa testowa",
        ]
        headerNamesEng = [
            "best score",
        ]
        headerNamesPl = [
            "Najlepszy wynik"
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()
        headerNames = headerNamesEng.copy() if not polish else headerNamesPl.copy()
        headerNames.extend(paramNames)
        # headerNames.append("visualization")
        headers = []
        if heuristicStep == 0:
            step = len(headerNames)
        else:
            step = heuristicStep


        for headerStart in range(0,len(headerNames),step):
            headerContent = headerPrefix.copy()
            headerContent.extend(headerNames[headerStart:headerStart+step])
            headerBegin = "\\begin{tabular}{|"+"|".join(["c" for _ in headerContent])+"|}\n"
            headerEnd = " \\\\\n"
            hline = "\\hline\n"
            header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
            headers.append(header)
        return headers

    def generateLatexGraphHeader(self, graphStep=2, polish=True):
        headers = []
        headerPrefixEng = [
            "test group",
            "graph group",
        ]
        headerPrefixPl = [
            "Grupa testowa",
            "Grupa sieci"
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()

        headerNames = GraphStats.vectorStatNamesPolish()
        if graphStep == 0:
            step = len(headerNames)
        else:
            step = graphStep


        for headerStart in range(0,len(headerNames),step):
            headerContent = headerPrefix.copy()
            headerContent.extend(headerNames[headerStart:headerStart+step])
    
            headerBegin = "\\begin{tabular}{|"+"|".join(["c" for header in headerContent])+"|}\n"
            headerEnd = " \\\\\n"
            hline = "\\hline\n"
            header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
            headers.append(header)
        return headers




    def generateLatexHeuristicRows(self, results, headers, heuristicStep=2):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        
        rowSegments = ["" for header in headers]
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        hline = "\\hline\n"
    
    
        for test_group in results:
            groupNames = test_group.split(" vs ")
    
            rowEnd = " \\\\\n"
            rowValueList = []
            rowPrefix = [
                    genCell(test_group, 1)
            ]
            rowValueList.append(genCell(results[test_group].best_value, 1))
            for param in paramNames:
                value = results[test_group].best_params[param] if param in results[test_group].best_params else "NA"
                rowValueList.append(genCell(value , 1))

            if heuristicStep == 0:
                step = len(rowValueList)
            else:
                step = heuristicStep


            for segmentIdx, rowStart in enumerate(range(0,len(rowValueList),step)):
                rowContent = rowPrefix.copy()
                rowContent.extend(rowValueList[rowStart:rowStart+step])
        
                row = " & ".join(rowContent) + rowEnd
                row += hline
                rowSegments[segmentIdx] = rowSegments[segmentIdx] + row
        

        return rowSegments

    def generateLatexGraphRows(self, results, headers, test_group_networks, distributionName, graphStep=2):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        rowSegments = ["" for header in headers]
        # print(f"len segments: {len(rowSegments)}")
    
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        genCline = lambda x,y: f"\\cline{{{x}-{y}}}\n"
        hline = "\\hline\n"
    
        for test_group in results:
            groupNames = test_group.split(" vs ")
    
            network1: ImmuneNetwork = test_group_networks[test_group][groupNames[0]]
            network2: ImmuneNetwork = test_group_networks[test_group][groupNames[1]]

            rowEnd = " \\\\\n"
            rowPrefix = [
                genCell(test_group, 2),
                genCell(groupNames[0], 1)
            ]
            rowValueList = []
            # row += genCell(img1, 1)
            for vectorStat in GraphStats(network1).toStatVector():
                rowValueList.append(genCell(vectorStat, 1))
    
            if graphStep == 0:
                step = len(rowValueList)
            else:
                step = graphStep
            # print(len(rowValueList))
            for segmentIdx, rowStart in enumerate(range(0,len(rowValueList),step)):
        
                rowContent = rowPrefix.copy()
                rowContent.extend(rowValueList[rowStart:rowStart+step])
        
                row = " & ".join(rowContent) + rowEnd
                row += genCline(2,len(rowContent))
                rowSegments[segmentIdx] = rowSegments[segmentIdx] + row
            rowValueList.clear()
    
            rowPrefix = [
                " ",
                genCell(groupNames[1], 1)
            ]
            # row += genCell(img2, 1)
            for vectorStat in GraphStats(network2).toStatVector():
                rowValueList.append(genCell(vectorStat, 1))
    
            for segmentIdx, rowStart in enumerate(range(0,len(rowValueList),step)):
                rowContent = rowPrefix.copy()
                rowContent.extend(rowValueList[rowStart:rowStart+step])
        
                row = " & ".join(rowContent) + rowEnd
                row += hline
                rowSegments[segmentIdx] = rowSegments[segmentIdx] + row
            rowValueList.clear()
        return rowSegments

    def printTable(self, result, test_group_networks, distributionName, actualDistributionName, heuristicStep=2, graphStep=2):
        best_params = [result[test_group].best_params for test_group in result]
        paramNames = set().union(*best_params)

        if heuristicStep == 0 and graphStep == 0:
            print("\\begin{landscape}")
        print(f"\\subsection{{Tabele dla {distributionName}}}")
        # print("\\vspace{1em}")

        sectionHeader = ""
        if heuristicStep == 0 and graphStep == 0:
            sectionHeader = "Tabela opisująca"
        else:
            sectionHeader = "Tabele opisujące"


        heuristicHeaders = self.generateLatexHeuristicHeader(result, heuristicStep)
        heuristicRowSegments = self.generateLatexHeuristicRows(result, heuristicHeaders, heuristicStep)


        print(f"\\subsubsection{{{sectionHeader} parametry tworzenia sieci znalezione przez heurystykę:}}")
        # print("\\vspace{1em}")


        for header, rowSegment in zip(heuristicHeaders, heuristicRowSegments):
            heuristicLatex = "\\begin{table}[htbp]\n"
            heuristicLatex += "\\begin{adjustbox}{max width =\\linewidth}\n"
            heuristicLatex = heuristicLatex + header
            heuristicLatex = heuristicLatex + rowSegment
            heuristicLatex += "\\end{tabular}\n"
            heuristicLatex += "\\end{adjustbox}\n"
            heuristicLatex += "\\end{table}\n"

            print()
            print("\\noindent")
            print(heuristicLatex)
            print()
            print("\\vspace{1em} \\\\")

        graphHeaders = self.generateLatexGraphHeader(graphStep)
        graphRowSegments = self.generateLatexGraphRows(result, graphHeaders, test_group_networks, distributionName, graphStep)

        print(f"\\subsubsection{{{sectionHeader} sieci wygenerowane na podstawie parametrów dla różnych grup}}")
        # print("\\vspace{1em}")

        for header, rowSegment in zip(graphHeaders, graphRowSegments):
            graphLatex = "\\begin{table}[htbp]\n"
            graphLatex += "\\begin{adjustbox}{max width =\\linewidth}\n"
            graphLatex = graphLatex + header
            graphLatex = graphLatex + rowSegment
            graphLatex += "\\end{tabular}\n"
            graphLatex += "\\end{adjustbox}\n"
            graphLatex += "\\end{table}\n"

            print()
            print("\\noindent")
            print(graphLatex)
            print("\\vspace{1em} \\\\")
        print()
        print("\\vspace{1em}")
        print()
        if heuristicStep == 0 and graphStep == 0:
            print("\\end{landscape}")

        print("\\subsubsection{Wizualizacje porównujące sieci wygenerowane dla różnych grup}")
        for test_group in result:
            figure = f'''
\\begin{{figure}}[H]
  \\centering
  \\includegraphics[width=1.0\\textwidth]{{figures/{actualDistributionName}_{test_group}.png}}
\\caption{{opis}}
  \\label{{fig:{actualDistributionName}_{test_group.replace(" ", "_")}}}
\\end{{figure}} \\\\
            '''
            print(figure)
            print()
            print("\\vspace{1em}")
            print()
            



