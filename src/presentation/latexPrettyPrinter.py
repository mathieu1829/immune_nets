from tqdm import tqdm

from src.entities import GraphStats, ImmuneNetwork
from src.creation.algorithms.simpleDistance import simpleDistance
from src.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from src.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from src.creation.distance.alignment import sequenceAligner
from src.creation.distance.levenshtein import levenshteinDistance 

class LatexPrettyPrinter:
    def generateLatexHeuristicHeader(self, results):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
        paramNames = [ name.replace("_", " ") for name in paramNames]
    
        header = ""
        headerNames = [
            "test group",
            "best score",
        ]
        headerNames.extend(paramNames)
        # headerNames.append("visualization")
    
        headerBegin = "\\begin{tabular}{|"+"|".join(["c" for header in headerNames])+"|}\n"
        headerEnd = " \\\\\n"
        hline = "\\hline\n"
        header = headerBegin + hline + " & ".join(headerNames) + headerEnd + hline
        return header

    def generateLatexGraphHeader(self, results):
        headers = []
        headerPrefix = [
            "test group",
            "graph group"
        ]
        headerNames = GraphStats.vectorStatNames()
        step = 4
        for headerStart in range(0,len(headerNames),step):
            header = ""
            headerContent = headerPrefix.copy()
            headerContent.extend(headerNames[headerStart:headerStart+step])
    
            headerBegin = "\\begin{tabular}{|"+"|".join(["c" for header in headerContent])+"|}\n"
            headerEnd = " \\\\\n"
            hline = "\\hline\n"
            header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
            headers.append(header)
        return headers




    def generateLatexHeuristicRows(self, results, allRepertoires):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        rows = ""
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        hline = "\\hline\n"
    
    
        for test_group in tqdm(results):
            groupNames = test_group.split(" vs ")
    
            # fig1, ax1 = plt.subplots(figsize=(6, 6))
            # fig2, ax2 = plt.subplots(figsize=(6, 6))
    
            # _ = graphVisualization(allRepertoires[test_group][groupNames[0]][0],network1, ax1)
            # _ = graphVisualization(allRepertoires[test_group][groupNames[0]][1],network2, ax2)
    
            # img1 = fig_to_base64(fig1)
            # img2 = fig_to_base64(fig2)
    
    
            rowEnd = " \\\\\n"
            rowValueList = []
            rowValueList.append(genCell(test_group, 1))
            rowValueList.append(genCell(results[test_group].best_value, 1))
            for param in paramNames:
                value = results[test_group].best_params[param] if param in results[test_group].best_params else "NA"
                rowValueList.append(genCell(value , 1))
        
            if rows == "":
                rows = " & ".join(rowValueList) + rowEnd
            else:
                rows += " & ".join(rowValueList) + rowEnd
                rows += hline

        return rows

    def generateLatexGraphRows(self, results, allRepertoires, headers, test_group_networks, genPlots=False):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        rowSegments = ["" for header in headers]
        # print(f"len segments: {len(rowSegments)}")
    
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        genCline = lambda x,y: f"\\cline{{{x}-{y}}}\n"
        hline = "\\hline\n"
    
        for test_group in tqdm(results):
            groupNames = test_group.split(" vs ")
    
            network1: ImmuneNetwork = test_group_networks[test_group][groupNames[0]]
            network2: ImmuneNetwork = test_group_networks[test_group][groupNames[1]]

            if genPlots:

    
            rowEnd = " \\\\\n"
            rowPrefix = [
                genCell(test_group, 2),
                genCell(groupNames[0], 1)
            ]
            rowValueList = []
            # row += genCell(img1, 1)
            for vectorStat in GraphStats(network1).toStatVector():
                rowValueList.append(genCell(vectorStat, 1))
    
            step = 4
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

    def printTable(self, result, allRepertoires, test_group_networks, genPlots=False):
        best_params = [result[test_group].best_params for test_group in result]
        paramNames = set().union(*best_params)

        heuristicLatex = ""
        heuristicLatex += self.generateLatexHeuristicHeader(result)
        heuristicLatex += self.generateLatexHeuristicRows(result, allRepertoires)
        heuristicLatex += "\\end{tabular}"

        print()
        print(heuristicLatex)
        print()

        graphHeaders = self.generateLatexGraphHeader(result)
        graphRowSegments = self.generateLatexGraphRows(result, allRepertoires, graphHeaders, test_group_networks, genPlots)

        for header, rowSegment in zip(graphHeaders, graphRowSegments):
            graphLatex = ""
            graphLatex = graphLatex + header
            graphLatex = graphLatex + rowSegment
            graphLatex += "\\end{tabular}"

            print()
            print(graphLatex)


