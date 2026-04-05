from immune_nets.entities import GraphStats, ImmuneNetwork
from immune_nets.creation.algorithms.simpleDistance import simpleDistance
from immune_nets.creation.algorithms.simpleBetaDistance import simpleBetaDistance
from immune_nets.creation.algorithms.simpleVectorBetaDistance import simpleVectorBetaDistance
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.distance.levenshtein import levenshteinDistance 
from immune_nets.analysis.visualization.multiGraphChart import multiGraphChart

class LatexPrettyPrinterTransposed:
    def generateLatexHeuristicHeader(self, results, polish=True):
        paramNames = [ test_group.replace("_","\\_") for test_group in results]
        
    
        headerPrefixEng = [
            "test group",
        ]
        headerPrefixPl = [
            "Grupa testowa",
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()
        headerNames = []
        headerNames.extend(paramNames)
        header = ""

        bold = lambda x: f"\\textbf{{{x}}}"


        headerContent = headerPrefix.copy()
        headerContent.extend(headerNames)
        headerContent = [ bold(header) for header in headerContent]
        headerBegin = "\\begin{tabular}{|"+"|".join(["c" for _ in headerContent])+"|}\n"
        headerEnd = " \\\\\n"
        hline = "\\hline\n"
        header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
        return header

    def generateLatexHeuristicRows(self, results, polish=True):

        best_params = [results[test_group].best_params for test_group in results]
        allParamNames = set().union(*best_params)
        paramNamesDict = {
            "algorithm_name":"Nazwa algorytmu",
            "distance":"Funckja dystansu",
            "substitution_matrix":"Macierz podstawień",
            "threshold":"Próg podobieństwa"
        }
        
        if polish:
            paramNames = [ paramNamesDict[param] for param in allParamNames]
        else:
            paramNames = [ name.replace("_", " ") for name in allParamNames]
    
        headerPrefixEng = [
        ]
        headerPrefixPl = [
        ]
        headerNamesEng = [
            "best score",
        ]
        headerNamesPl = [
            "Najlepszy wynik"
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()
        headerNamesInitial = headerNamesEng.copy() if not polish else headerNamesPl.copy()
        headerNames = headerPrefix.copy()
        headerNames.extend(headerNamesInitial)
        headerNames.extend(paramNames)

        rows = ""
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        hline = "\\hline\n"
        rowEnd = " \\\\\n"
        bold = lambda x: f"\\textbf{{{x}}}"

        rowData = [ [genCell(bold(header), 1)] for header in headerNames]

        for test_group in results:
            result = results[test_group]
            rowData[0].append(genCell(result.best_value,1))

            for paramIdx, paramName in enumerate(allParamNames):
                value = results[test_group].best_params[paramName] if paramName in results[test_group].best_params else "NA"
                rowData[paramIdx+1].append(genCell(value, 1))

        for rowContent in rowData:
            row = " & ".join(rowContent) + rowEnd
            row += hline
            rows = rows + row
        

        return rows

    def generateLatexGraphHeader(self, results, polish=True):
        multicolumn = lambda content, x, firstColumn=False: f"\\multicolumn{{{x}}}{{{'|c|' if firstColumn else 'c|'}}}{{{content}}}"
        paramNames = [ test_group.replace("_","\\_") for test_group in results]
        headers = []
        headerPrefixEng = [
            "test group",
        ]
        headerPrefixPl = [
            "Grupa testowa",
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()

        headerNames = []
        headerNames.extend(paramNames)

        bold = lambda x: f"\\textbf{{{x}}}"

        headerContent = headerPrefix.copy()
        headerContent.extend(headerNames)
        headerContent = [multicolumn(bold(header), 2) if headerIdx != 0 else multicolumn(bold(header), 1, firstColumn=True) for headerIdx, header in enumerate(headerContent)]

        headerBegin = "\\begin{tabular}{|"+"|".join(["c" for header in range(2*len(GraphStats.vectorStatNamesPolish())+1)])+"|}\n"
        headerEnd = " \\\\\n"
        hline = "\\hline\n"
        header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
        return header

    def generateLatexGraphRows(self, results, test_group_networks, distributionName, polish=True):
        best_params = [results[test_group].best_params for test_group in results]
        paramNames = set().union(*best_params)
    
        # print(f"len segments: {len(rowSegments)}")
    
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        genCline = lambda x,y: f"\\cline{{{x}-{y}}}\n"
        hline = "\\hline\n"
        rowEnd = " \\\\\n"
        bold = lambda x: f"\\textbf{{{x}}}"

        headerPrefixEng = [
            "graph group",
        ]
        headerPrefixPl = [
            "Grupa sieci"
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()

        headerNames = headerPrefix.copy()
        headerNames.extend(GraphStats.vectorStatNamesPolish())

        rowData = [ [bold(headerName)] for headerName in headerNames]

        for test_group in results:
            result = results[test_group]
            groupNames = test_group.split(" vs ")
    
            network1: ImmuneNetwork = test_group_networks[test_group][groupNames[0]]
            network2: ImmuneNetwork = test_group_networks[test_group][groupNames[1]]

            stats1 = GraphStats(network1).toStatVector()
            stats2 = GraphStats(network2).toStatVector()

            for idx, (stat1, stat2) in enumerate(zip( stats1, stats2)):
                if idx == 0:
                    rowData[0].append(genCell(groupNames[0],1))
                    rowData[0].append(genCell(groupNames[1],1))
                rowData[idx+1].append(genCell(stat1,1))
                rowData[idx+1].append(genCell(stat2,1))


        rows = ""
        for rowRecord in rowData:
            row = " & ".join(rowRecord) + rowEnd
            row += hline
            rows = rows + row

        return rows

       
    


    def printTable(self, result, test_group_networks, distributionName, actualDistributionName):

        # print("\\vspace{1em}")

        sectionHeader = "Tabela opisująca"


        heuristicHeaders = self.generateLatexHeuristicHeader(result)
        heuristicRows = self.generateLatexHeuristicRows(result)


        heuristicLatex = "\\begin{table}[htbp]\n"
        heuristicLatex += f"\\caption{{{sectionHeader} parametry tworzenia sieci znalezione przez heurystykę dla {distributionName}}}\n"
        heuristicLatex += "\\begin{adjustbox}{max width =\\linewidth}\n"
        heuristicLatex = heuristicLatex + heuristicHeaders
        heuristicLatex = heuristicLatex + heuristicRows
        heuristicLatex += "\\end{tabular}\n"
        heuristicLatex += "\\end{adjustbox}\n"
        heuristicLatex += "\\end{table}\n"

        print()
        print("\\noindent")
        print(heuristicLatex)
        print()
        print("\\vspace{1em}")

        graphHeader = self.generateLatexGraphHeader(result)
        graphRows = self.generateLatexGraphRows(result, test_group_networks, distributionName)

        # print("\\vspace{1em}")

        graphLatex = "\\begin{table}[htbp]\n"
        graphLatex += f"\\caption{{{sectionHeader} sieci dla różnych grup wygenerowane na podstawie parametrów zopytmalizowanych dla {distributionName}}}\n"
        graphLatex += "\\begin{adjustbox}{max width =\\linewidth}\n"
        graphLatex = graphLatex + graphHeader
        graphLatex = graphLatex + graphRows
        graphLatex += "\\end{tabular}\n"
        graphLatex += "\\end{adjustbox}\n"
        graphLatex += "\\end{table}\n"

        print()
        print("\\noindent")
        print(graphLatex)
        print("\\vspace{1em}")

        print()
        print("\\vspace{1em}")
        print()

        for test_group in result:
            groupNames = test_group.split(" vs ")
            figure = f'''
\\begin{{figure}}[H]
  \\centering
  \\includegraphics[width=1.0\\textwidth]{{figures/{actualDistributionName}_{test_group}.png}}
\\caption{{Wizualizacja porównująca przykładowe sieci dla grup {groupNames[0].replace("_", "\\_")} i {groupNames[1].replace("_", "\\_")} dla {distributionName}. }}
  \\label{{fig:{actualDistributionName}_{test_group.replace(" ", "_")}}}
\\end{{figure}}
            '''
            print(figure)
            print()
            print("\\vspace{1em}")
            print()
            
    def generateLatexOptHeuristicHeader(self, results, polish=True):
        
    
        headerPrefixEng = [
            "Heuristic run",
        ]
        headerPrefixPl = [
            "Numer uruchomienia heurystyki",
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()
        header = ""

        bold = lambda x: f"\\textbf{{{x}}}"


        headerContent = headerPrefix.copy()
        headerContent.extend([ str(i+1) for i in range(len(results))])
        headerContent = [ bold(header) for header in headerContent]
        headerBegin = "\\begin{tabular}{|"+"|".join(["c" for _ in headerContent])+"|}\n"
        headerEnd = " \\\\\n"
        hline = "\\hline\n"
        header = headerBegin + hline + " & ".join(headerContent) + headerEnd + hline
        return header

    def generateLatexOptHeuristicRows(self, results, polish=True):

        best_params = [result.params for result in results]
        allParamNames = set().union(*best_params)
        user_attrs = [result.user_attrs for result in results]
        allUserAttrs = set().union(*user_attrs)
        paramNamesDict = {
            "algorithm_name":"Nazwa algorytmu",
            "distance":"Funckja dystansu",
            "substitution_matrix":"Macierz podstawień",
            "threshold":"Próg podobieństwa"
        }
        userAttrsDict = {
            'healthy vs covid score': "Średni dystans dla grup healthy vs covid",
            'covid_1 vs covid_2 score': "Średni dystans dla grup covid_1 vs covid_2",
            'healthy_1 vs healthy_2 score': "Średni dystans dla grup healthy_1 vs healthy_2",
            'between': "Średni dystans między grupami",
            'within': "Średni dystans wewnątrz grup"
        }
        
        if polish:
            paramNames = [ paramNamesDict[param] for param in allParamNames]
            userAttrsNames = [ userAttrsDict[attrName] for attrName in allUserAttrs]
        else:
            paramNames = [ name.replace("_", " ") for name in allParamNames]
            userAttrsNames = [ attrName.replace("_", " ") for attrName in allUserAttrs]
    
        headerPrefixEng = [
        ]
        headerPrefixPl = [
        ]
        headerNamesEng = [
            "best score",
        ]
        headerNamesPl = [
            "Najlepszy wynik"
        ]
        headerPrefix = headerPrefixEng.copy() if not polish else headerPrefixPl.copy()
        headerNamesInitial = headerNamesEng.copy() if not polish else headerNamesPl.copy()
        headerNames = headerPrefix.copy()
        headerNames.extend(headerNamesInitial)
        headerNames.extend(userAttrsNames)
        headerNames.extend(paramNames)

        rows = ""
        formatValue = lambda value: f"{value:.2f}" if isinstance(value, float) else f"{value.replace("_", "\\_")}"
        genCell = lambda value,x: f"{formatValue(value)}" if x == 1 else "\\multirow{"+f"{x}"+"}{*}{"+f"{formatValue(value)}"+"}"
        hline = "\\hline\n"
        rowEnd = " \\\\\n"
        bold = lambda x: f"\\textbf{{{x}}}"

        rowData = [ [genCell(bold(header), 1)] for header in headerNames]

        for result in results:
            
            rowData[0].append(genCell(result.values[0],1))

            for attrIdx, attrName in enumerate(allUserAttrs):
                value = result.user_attrs[attrName]
                rowData[attrIdx+1].append(genCell(value, 1))

            for paramIdx, paramName in enumerate(allParamNames):
                value = result.params[paramName] if paramName in result.params else "NA"
                rowData[paramIdx+len(result.user_attrs)+1].append(genCell(value, 1))

        for rowContent in rowData:
            row = " & ".join(rowContent) + rowEnd
            row += hline
            rows = rows + row


        return rows

    def printOptTable(self, result,  distributionName):

        # print("\\vspace{1em}")

        sectionHeader = "Tabela opisująca"


        heuristicHeaders = self.generateLatexOptHeuristicHeader(result)
        heuristicRows = self.generateLatexOptHeuristicRows(result)

        heuristicLatex = "\\begin{table}[htbp]\n"
        heuristicLatex += f"\\caption{{{sectionHeader} parametry tworzenia sieci znalezione przez heurystykę dla {distributionName}}}\n"
        heuristicLatex += "\\begin{adjustbox}{max width =\\linewidth}\n"
        heuristicLatex = heuristicLatex + heuristicHeaders
        heuristicLatex = heuristicLatex + heuristicRows
        heuristicLatex += "\\end{tabular}\n"
        heuristicLatex += "\\end{adjustbox}\n"
        heuristicLatex += "\\end{table}\n"

        print()
        print("\\noindent")
        print(heuristicLatex)
        print()
        print("\\vspace{1em}")

