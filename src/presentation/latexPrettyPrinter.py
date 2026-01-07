from src.entities import GraphStats, ImmuneNetwork
from src.creation.algorithms.simpleDistance import simpleDistance
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 
from src.analysis.visualization.multiGraphChart import multiGraphChart

class LatexPrettyPrinter:
    def generateLatexHeuristicHeader(self, results, polish=True):
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
        step = 3 

        for headerStart in range(0,len(headerNames),step):
            headerContent = headerPrefix.copy()
            headerContent.extend(headerNames[headerStart:headerStart+step])

            headerBegin = "\\begin{tabular}{|"+"|".join(["c" for _ in headerContent])+"|}\n"
            headerEnd = " \\\\\n"
            hline = "\\hline\n"
            header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
            headers.append(header)
        return headers

    def generateLatexGraphHeader(self, polish=True):
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
        step = 3
        for headerStart in range(0,len(headerNames),step):
            headerContent = headerPrefix.copy()
            headerContent.extend(headerNames[headerStart:headerStart+step])
    
            headerBegin = "\\begin{tabular}{|"+"|".join(["c" for header in headerContent])+"|}\n"
            headerEnd = " \\\\\n"
            hline = "\\hline\n"
            header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
            headers.append(header)
        return headers




    def generateLatexHeuristicRows(self, results, headers):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        
        rowSegments = ["" for header in headers]
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value}"
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

            step = 3
            for segmentIdx, rowStart in enumerate(range(0,len(rowValueList),step)):
                rowContent = rowPrefix.copy()
                rowContent.extend(rowValueList[rowStart:rowStart+step])
        
                row = " & ".join(rowContent) + rowEnd
                row += hline
                rowSegments[segmentIdx] = rowSegments[segmentIdx] + row
        

        return rowSegments

    def generateLatexGraphRows(self, results, headers, test_group_networks, distributionName, genPlots=False):
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

            # if genPlots:
            #     plotTitles = [ f"Sieć próbki {group} dla {distributionName}" for group in groupNames] 
            #     multiGraphChart(plotTitles, [network1, network2], f"{distributionName}_{test_group}.png")

    
            rowEnd = " \\\\\n"
            rowPrefix = [
                genCell(test_group, 2),
                genCell(groupNames[0], 1)
            ]
            rowValueList = []
            # row += genCell(img1, 1)
            for vectorStat in GraphStats(network1).toStatVector():
                rowValueList.append(genCell(vectorStat, 1))
    
            step = 3
            # print(len(rowValueList))
            for segmentIdx, rowStart in enumerate(range(0,len(rowValueList),step)):
                # print(f"rowStart: {rowStart}")
                # print(f"segmentIdx: {segmentIdx}")
        
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
                # row += genCline(1,len(rowContent)+1)
                row += hline
                rowSegments[segmentIdx] = rowSegments[segmentIdx] + row
            rowValueList.clear()
            # rows += " & ".join(rowValueList) + rowEnd
            # rows += genCline(1,len(rowValueList)+1)
        return rowSegments

    def printTable(self, result, test_group_networks, distributionName, genPlots=False):
        best_params = [result[test_group].best_params for test_group in result]
        paramNames = set().union(*best_params)


        print(f"Tabele dla {distributionName}")
        print("\\vspace{1em}")

        heuristicHeaders = self.generateLatexHeuristicHeader(result)
        heuristicRowSegments = self.generateLatexHeuristicRows(result, heuristicHeaders)


        print(f"Tabele opisujące parametry tworzenia sieci znalezione przez heurystykę")
        print("\\vspace{1em}")


        for header, rowSegment in zip(heuristicHeaders, heuristicRowSegments):
            heuristicLatex = ""
            heuristicLatex = heuristicLatex + header
            heuristicLatex = heuristicLatex + rowSegment
            heuristicLatex += "\\end{tabular}"

            print()
            print("\\noindent")
            print(heuristicLatex)
            print()
            print("\\vspace{1em}")

        graphHeaders = self.generateLatexGraphHeader(result)
        graphRowSegments = self.generateLatexGraphRows(result, graphHeaders, test_group_networks, distributionName, genPlots)

        print("Tabele porównujące sieci wygenerowane na podstawie parametrów dla różnych grup")
        print("\\vspace{1em}")

        for header, rowSegment in zip(graphHeaders, graphRowSegments):
            graphLatex = ""
            graphLatex = graphLatex + header
            graphLatex = graphLatex + rowSegment
            graphLatex += "\\end{tabular}"

            print()
            print("\\noindent")
            print(graphLatex)
            print("\\vspace{1em}")
        print()
        print("\\vspace{1em}")
        print()



