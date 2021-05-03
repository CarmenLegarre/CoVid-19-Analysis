#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 09:36:15 2021

@author: carmen
"""

import wget
import os
import zipfile

def download_file(tag, url):
    """
    The function download data regarding the given tag, and saves the 
    information in a folder with a name equal to the tag
    """
    try:
        wget.download(url, '/home/carmen/Repositories/CoVid/CoVid-19-Analysis/code/Covid_data_Spain/data/'+tag) 
    except:
        print("Impossible to download data about " +tag)
        
def remove_create(tag):
    """
    It removes data from the given tag folder to avoid duplicated information.
    """
    os.chdir(tag)
    for document in os.listdir("./"):
        os.remove(document)
    os.chdir("..")
    
    

    
# Comunidad autonoma - abreviatura ISO
    
# Andalucia - AN
tag = "AN"
remove_create(tag)
url = 'https://www.juntadeandalucia.es/institutodeestadisticaycartografia/badea/stpivot/stpivot/Print?cube=0b7991ef-de75-4e76-bf77-1a53e019ab33&type=0&foto=si&ejecutaDesde=&codConsulta=39409&consTipoVisua=JP'
download_file(tag, url)

# Aragon - AR
tag = "AR"
remove_create(tag)
url = 'https://www.aragon.es/documents/20127/38742837/casos_coronavirus_aragon.csv'
download_file(tag, url)

# Asturias - AS
tag = "AS"
remove_create(tag)
# este url no es correcto. Asturias parece no tener la posibilidad de descargar 
# datos
url = 'https://www.asturias.es/portal/site/webasturias/template.PAGE/BuscadorOpenDataTransparencia/?accion=resultados_view&busqueda=&busq-tipo=7b3106c9a1750610VgnVCM10000077020a0aRCRD&FechaCheckHidden1=F-FIN&FechaCheckCompHidden1=1&canales=&fechaInicio=&fechaFin=&vgnextoid=320afaf18a869510VgnVCM100000ce212b0aRCRD&i18n.http.lang=es'
download_file(tag, url)

# Carnarias - CN
# fuente no oficial
tag = "CN"
remove_create(tag)
url = 'data:application/octet-stream;charset=utf-8,%EF%BB%BFedad%2CMujer%2CHombre%2CTotal%2C%25%0A0-9%20a%C3%B1os%2C566%2C596%2C1162%2C5.18%0A10-19%20a%C3%B1os%2C931%2C1048%2C1979%2C8.82%0A20-29%20a%C3%B1os%2C1879%2C1896%2C3775%2C16.82%0A30-39%20a%C3%B1os%2C1924%2C1822%2C3746%2C16.69%0A40-49%20a%C3%B1os%2C2027%2C1939%2C3966%2C17.67%0A50-59%20a%C3%B1os%2C1666%2C1603%2C3269%2C14.56%0A60-69%20a%C3%B1os%2C1075%2C1077%2C2152%2C9.59%0A70-79%20a%C3%B1os%2C672%2C718%2C1390%2C6.19%0A80-89%20a%C3%B1os%2C455%2C349%2C804%2C3.58%0A%3E%3D90%20a%C3%B1os%2C147%2C56%2C203%2C0.9%0ADesconocida%2C0%2C0%2C0%2C0'
download_file(tag, url)

# Cantabria - CB
tag = "CB"
remove_create(tag)
url = 'https://serviweb.scsalud.es:10443/ficheros/COVID19_historico.csv'
download_file(tag, url)

# Castilla la Mancha - CM
tag = "CM"
remove_create(tag)
# informacion de numeros de PCR
url = 'data:application/octet-stream;charset=utf-8,%EF%BB%BFD%C3%ADa%20de%20la%20semana%2CCasos%20PCR%0A01%2F05%2C57%0A02%2F05%2C50%0A03%2F05%2C33%0A04%2F05%2C30%0A05%2F05%2C40%0A06%2F05%2C64%0A07%2F05%2C40%0A08%2F05%2C53%0A09%2F05%2C41%0A10%2F05%2C65%0A11%2F05%2C44%0A12%2F05%2C44%0A13%2F05%2C22%0A14%2F05%2C17%0A15%2F05%2C43%0A16%2F05%2C74%0A17%2F05%2C31%0A18%2F05%2C30%0A19%2F05%2C29%0A20%2F05%2C62%0A21%2F05%2C50%0A22%2F05%2C41%0A23%2F05%2C25%0A24%2F05%2C34%0A25%2F05%2C20%0A26%2F05%2C43%0A27%2F05%2C30%0A28%2F05%2C84%0A29%2F05%2C54%0A30%2F05%2C62%0A31%2F05%2C50%0A01%2F06%2C25%0A02%2F06%2C51%0A03%2F06%2C45%0A04%2F06%2C38%0A05%2F06%2C38%0A06%2F06%2C41%0A07%2F06%2C24%0A08%2F06%2C37%0A09%2F06%2C68%0A10%2F06%2C39%0A11%2F06%2C62%0A12%2F06%2C38%0A13%2F06%2C36%0A14%2F06%2C27%0A15%2F06%2C8%0A16%2F06%2C17%0A17%2F06%2C17%0A18%2F06%2C49%0A19%2F06%2C17%0A20%2F06%2C29%0A21%2F06%2C25%0A22%2F06%2C22%0A23%2F06%2C17%0A24%2F06%2C26%0A25%2F06%2C24%0A26%2F06%2C36%0A27%2F06%2C21%0A28%2F06%2C16%0A29%2F06%2C9%0A30%2F06%2C30%0A1%2F07%2C23%0A2%2F07%2C18%0A3%2F07%2C35%0A6%2F07%2C37%0A7%2F07%2C22%0A8%2F07%2C14%0A9%2F07%2C30%0A10%2F07%2C20%0A13%2F07%2C23%0A14%2F07%2C18%0A15%2F07%2C10%0A16%2F07%2C40%0A17%2F07%2C36%0A20%2F07%2C52%0A21%2F07%2C14%0A22%2F07%2C25%0A23%2F07%2C27%0A29%2F07%2C32%0A30%2F07%2C56%0A31%2F07%2C68%0A03%2F08%2C146%0A04%2F08%2C50%0A05%2F08%2C81%0A06%2F08%2C85%0A07%2F08%2C120%0A10%2F08%2C384%0A11%2F08%2C165%0A12%2F08%2C180%0A13%2F08%2C186%0A14%2F08%2C201%0A17%2F08%2C548%0A18%2F08%2C138%0A19%2F08%2C313%0A20%2F08%2C343%0A21%2F08%2C370%0A24%2F08%2C893%0A25%2F08%2C408%0A26%2F08%2C387%0A27%2F08%2C486%0A28%2F08%2C507%0A31%2F08%2C1423%0A01%2F09%2C490%0A02%2F09%2C616%0A03%2F09%2C628%0A04%2F09%2C769%0A07%2F09%2C1845%0A08%2F09%2C541%0A09%2F09%2C518%0A10%2F09%2C798%0A11%2F09%2C686%0A14%2F09%2C2168%0A15%2F09%2C463%0A16%2F09%2C700%0A17%2F09%2C817%0A18%2F09%2C978%0A21%2F09%2C2470%0A22%2F09%2C422%0A23%2F09%2C747%0A24%2F09%2C887%0A25%2F09%2C891%0A28%2F09%2C1940%0A29%2F09%2C384%0A30%2F09%2C566%0A01%2F10%2C828%0A02%2F10%2C642%0A05%2F10%2C1608%0A06%2F10%2C402%0A07%2F10%2C670%0A08%2F10%2C649%0A09%2F10%2C695%0A13%2F10%2C2120%0A14%2F10%2C365%0A15%2F10%2C746%0A16%2F10%2C909%0A19%2F10%2C2026%0A20%2F10%2C537%0A21%2F10%2C991%0A22%2F10%2C1031%0A23%2F10%2C982%0A26%2F10%2C2773%0A27%2F10%2C621%0A28%2F10%2C1212%0A29%2F10%2C1110%0A30%2F10%2C1138%0A31%2F10%2C1033%0A01%2F11%2C864%0A02%2F11%2C544%0A03%2F11%2C767%0A04%2F11%2C1041%0A05%2F11%2C1098%0A06%2F11%2C930%0A07%2F11%2C942%0A08%2F11%2C767%0A09%2F11%2C277%0A10%2F11%2C605%0A11%2F11%2C989%0A12%2F11%2C960%0A13%2F11%2C824%0A14%2F11%2C915%0A15%2F11%2C706%0A16%2F11%2C240%0A17%2F11%2C503%0A18%2F11%2C713%0A19%2F11%2C784%0A20%2F11%2C588%0A21%2F11%2C670%0A22%2F11%2C392%0A23%2F11%2C145%0A24%2F11%2C358%0A25%2F11%2C560%0A26%2F11%2C682%0A27%2F11%2C492%0A28%2F11%2C305%0A29%2F11%2C285%0A30%2F11%2C111%0A01%2F12%2C309%0A02%2F12%2C588%0A03%2F12%2C451%0A04%2F12%2C428%0A07%2F12%2C1123%0A09%2F12%2C714%0A10%2F12%2C309%0A11%2F12%2C474%0A14%2F12%2C1315%0A15%2F12%2C348%0A16%2F12%2C680%0A17%2F12%2C789%0A18%2F12%2C535%0A19%2F12%2C589%0A20%2F12%2C545%0A21%2F12%2C150%0A22%2F12%2C375%0A23%2F12%2C866%0A24%2F12%2C603%0A25%2F12%2C833%0A26%2F12%2C253%0A27%2F12%2C570%0A29%2F12%2C470%0A30%2F12%2C1132%0A31%2F12%2C958%0A01%2F01%2C539%0A02%2F01%2C601%0A03%2F01%2C498'
download_file(tag, url)
# informacion de hospitalizaciones
url = 'data:application/octet-stream;charset=utf-8,%EF%BB%BFD%C3%ADa%20de%20la%20semana%2CUCI%2CHospitalizaciones%0A28%2F03%2C289%2C2977%0A29%2F03%2C299%2C3018%0A30%2F03%2C302%2C3134%0A31%2F03%2C344%2C3225%0A1%2F04%2C353%2C3230%0A2%2F04%2C355%2C3184%0A3%2F04%2C355%2C3168%0A4%2F04%2C360%2C3133%0A5%2F04%2C357%2C2950%0A6%2F04%2C354%2C2901%0A7%2F04%2C360%2C2908%0A8%2F04%2C354%2C2731%0A9%2F04%2C342%2C2575%0A10%2F04%2C329%2C2393%0A11%2F04%2C316%2C2198%0A12%2F04%2C302%2C2067%0A13%2F04%2C300%2C2047%0A14%2F04%2C302%2C1973%0A15%2F04%2C294%2C1867%0A16%2F04%2C281%2C1799%0A17%2F04%2C271%2C1636%0A18%2F04%2C248%2C1625%0A19%2F04%2C239%2C1458%0A20%2F04%2C241%2C1430%0A21%2F04%2C230%2C1422%0A22%2F04%2C222%2C1336%0A23%2F04%2C214%2C1261%0A24%2F04%2C208%2C1172%0A25%2F04%2C202%2C1111%0A26%2F04%2C190%2C1026%0A27%2F04%2C196%2C989%0A28%2F04%2C189%2C995%0A29%2F04%2C181%2C950%0A30%2F04%2C180%2C878%0A01%2F05%2C169%2C771%0A02%2F05%2C165%2C713%0A03%2F05%2C150%2C713%0A04%2F05%2C142%2C694%0A05%2F05%2C140%2C688%0A06%2F05%2C132%2C679%0A07%2F05%2C122%2C671%0A08%2F05%2C118%2C619%0A09%2F05%2C114%2C588%0A10%2F05%2C105%2C513%0A11%2F05%2C99%2C523%0A12%2F05%2C98%2C505%0A13%2F05%2C93%2C460%0A14%2F05%2C91%2C444%0A15%2F05%2C90%2C394%0A16%2F05%2C83%2C370%0A17%2F05%2C80%2C327%0A18%2F05%2C80%2C323%0A19%2F05%2C74%2C303%0A20%2F05%2C69%2C265%0A21%2F05%2C62%2C226%0A22%2F05%2C59%2C206%0A23%2F05%2C50%2C180%0A24%2F05%2C50%2C181%0A25%2F05%2C47%2C187%0A26%2F05%2C44%2C177%0A27%2F05%2C45%2C165%0A28%2F05%2C41%2C165%0A29%2F05%2C36%2C153%0A30%2F05%2C33%2C142%0A31%2F05%2C34%2C138%0A01%2F06%2C33%2C130%0A02%2F06%2C32%2C131%0A03%2F06%2C30%2C113%0A04%2F06%2C30%2C115%0A05%2F06%2C31%2C102%0A06%2F06%2C29%2C85%0A07%2F06%2C28%2C85%0A08%2F06%2C29%2C84%0A09%2F06%2C28%2C82%0A10%2F06%2C28%2C66%0A11%2F06%2C25%2C56%0A12%2F06%2C25%2C56%0A13%2F06%2C23%2C54%0A14%2F06%2C21%2C52%0A15%2F06%2C21%2C55%0A16%2F06%2C20%2C51%0A17%2F06%2C19%2C53%0A18%2F06%2C19%2C52%0A19%2F06%2C17%2C55%0A20%2F06%2C17%2C46%0A21%2F06%2C16%2C48%0A25%2F06%2C15%2C44%0A26%2F06%2C12%2C43%0A27%2F06%2C13%2C46%0A29%2F06%2C12%2C48%0A30%2F06%2C11%2C51%0A1%2F07%2C10%2C46%0A2%2F07%2C9%2C45%0A3%2F07%2C8%2C44%0A6%2F07%2C7%2C57%0A7%2F07%2C7%2C52%0A8%2F07%2C9%2C48%0A9%2F07%2C8%2C41%0A13%2F07%2C8%2C36%0A14%2F07%2C8%2C43%0A15%2F07%2C7%2C44%0A16%2F07%2C9%2C41%0A17%2F07%2C10%2C44%0A20%2F07%2C11%2C39%0A21%2F07%2C10%2C39%0A22%2F07%2C9%2C45%0A23%2F07%2C9%2C42%0A30%2F07%2C9%2C35%0A31%2F07%2C8%2C32%0A03%2F08%2C10%2C51%0A04%2F08%2C10%2C53%0A05%2F08%2C10%2C45%0A06%2F08%2C9%2C46%0A07%2F08%2C9%2C46%0A10%2F08%2C10%2C44%0A11%2F08%2C9%2C56%0A12%2F08%2C9%2C67%0A13%2F08%2C8%2C67%0A14%2F08%2C8%2C65%0A17%2F08%2C8%2C84%0A18%2F08%2C10%2C98%0A19%2F08%2C10%2C92%0A20%2F08%2C10%2C102%0A21%2F08%2C9%2C119%0A24%2F08%2C10%2C157%0A25%2F08%2C17%2C166%0A26%2F08%2C20%2C175%0A27%2F08%2C22%2C180%0A28%2F08%2C23%2C186%0A31%2F08%2C25%2C251%0A1%2F09%2C29%2C256%0A2%2F09%2C33%2C261%0A3%2F09%2C39%2C277%0A4%2F09%2C36%2C311%0A7%2F09%2C40%2C358%0A8%2F09%2C38%2C372%0A9%2F09%2C36%2C390%0A10%2F09%2C36%2C408%0A11%2F09%2C41%2C419%0A14%2F09%2C50%2C479%0A15%2F09%2C50%2C492%0A16%2F09%2C50%2C511%0A17%2F09%2C60%2C519%0A18%2F09%2C63%2C506%0A21%2F09%2C66%2C547%0A22%2F09%2C68%2C554%0A23%2F09%2C72%2C549%0A24%2F09%2C70%2C557%0A25%2F09%2C70%2C575%0A28%2F09%2C73%2C577%0A29%2F09%2C66%2C539%0A30%2F09%2C69%2C536%0A01%2F10%2C74%2C529%0A02%2F10%2C74%2C508%0A05%2F10%2C77%2C554%0A06%2F10%2C76%2C534%0A07%2F10%2C78%2C506%0A08%2F10%2C77%2C481%0A09%2F10%2C73%2C489%0A13%2F10%2C75%2C561%0A14%2F10%2C74%2C550%0A15%2F10%2C78%2C573%0A16%2F10%2C75%2C560%0A19%2F10%2C71%2C608%0A20%2F10%2C74%2C597%0A21%2F10%2C70%2C608%0A22%2F10%2C74%2C636%0A23%2F10%2C72%2C653%0A26%2F10%2C76%2C693%0A27%2F10%2C76%2C685%0A28%2F10%2C79%2C714%0A29%2F10%2C77%2C746%0A30%2F10%2C86%2C771%0A31%2F10%2C83%2C774%0A01%2F11%2C89%2C848%0A02%2F11%2C87%2C871%0A03%2F11%2C81%2C838%0A04%2F11%2C82%2C850%0A05%2F11%2C84%2C837%0A06%2F11%2C84%2C803%0A07%2F11%2C87%2C795%0A08%2F11%2C91%2C804%0A09%2F11%2C87%2C840%0A10%2F11%2C96%2C817%0A11%2F11%2C101%2C792%0A12%2F11%2C103%2C772%0A13%2F11%2C106%2C772%0A14%2F11%2C105%2C698%0A15%2F11%2C104%2C732%0A16%2F11%2C101%2C778%0A17%2F11%2C97%2C756%0A18%2F11%2C97%2C714%0A19%2F11%2C96%2C647%0A20%2F11%2C101%2C635%0A21%2F11%2C97%2C591%0A22%2F11%2C97%2C562%0A23%2F11%2C93%2C573%0A24%2F11%2C91%2C542%0A25%2F11%2C90%2C534%0A26%2F11%2C94%2C503%0A27%2F11%2C92%2C498%0A28%2F11%2C87%2C506%0A29%2F11%2C88%2C519%0A30%2F11%2C94%2C537%0A01%2F12%2C95%2C514%0A02%2F12%2C96%2C499%0A03%2F12%2C93%2C506%0A04%2F12%2C93%2C478%0A07%2F12%2C89%2C463%0A09%2F12%2C83%2C449%0A10%2F12%2C76%2C475%0A11%2F12%2C76%2C476%0A14%2F12%2C70%2C510%0A15%2F12%2C67%2C487%0A16%2F12%2C63%2C483%0A17%2F12%2C65%2C466%0A18%2F12%2C66%2C446%0A21%2F12%2C61%2C451%0A22%2F12%2C63%2C462%0A23%2F12%2C67%2C447%0A28%2F12%2C82%2C466%0A29%2F12%2C85%2C505%0A04%2F01%2C89%2C597' 
download_file(tag, url)

# Castilla y Leo - CL
tag = "CL"
remove_create(tag)
# informacion sobre el numero de test realizados
url = 'https://analisis.datosabiertos.jcyl.es/explore/dataset/pruebas-realizados-coronavirus/download/?format=csv&disjunctive.provincia=true&refine.provincia=Castilla+y+Le%C3%B3n&timezone=Europe/Madrid&lang=es&use_labels_for_header=true&csv_separator=%3B'
download_file(tag, url)
# situacion epidemiologica
url = 'https://analisis.datosabiertos.jcyl.es/explore/dataset/situacion-epidemiologica-coronavirus-en-castilla-y-leon/download/?format=csv&timezone=Europe/Madrid&lang=es&use_labels_for_header=true&csv_separator=%3B'
download_file(tag, url)
# hospitalizaciones
url = 'https://analisis.datosabiertos.jcyl.es/explore/dataset/situacion-de-hospitalizados-por-coronavirus-en-castilla-y-leon/download/?format=csv&timezone=Europe/Madrid&lang=es&use_labels_for_header=true&csv_separator=%3B'
download_file(tag, url)

# Cataluña - CT
tag = "CT"
remove_create(tag)
url = 'https://analisi.transparenciacatalunya.cat/api/views/c7sd-zy9j/rows.csv?accessType=DOWNLOAD&sorting=true'
download_file(tag, url)

# Extremadura - EX
tag = "EX"
remove_create(tag)
url = ''
download_file(tag, url)

# Galicia - GA
tag = "GA"
remove_create(tag)
url = ''
download_file(tag, url)

# Islas Baleares - IB
tag = "IB"
url = ''
download_file(tag, url)

# La Rioja - RI
tag = "RI"
remove_create(tag)
# test
url = 'https://ias1.larioja.org/opendata/download?r=Y2Q9ODU2fGNmPTAz'
download_file(tag, url)
# positivos
url = 'https://ias1.larioja.org/opendata/download?r=Y2Q9ODU3fGNmPTAz'
download_file(tag, url)

# Comunidad de Madrid - MD
tag = "MD"
remove_create(tag)
url = ''
download_file(tag, url)

# Region de Murcia - MC
tag = "MC"
remove_create(tag)
url = 'https://datos.comunidad.madrid/catalogo/dataset/7da43feb-8d4d-47e0-abd5-3d022d29d09e/resource/b2a3a3f9-1f82-42c2-89c7-cbd3ef801412/download/covid19_tia_muni_y_distritos.csv'
download_file(tag, url)

# Comunidad Foral de Navarra - NC
tag = "NC"
remove_create(tag)
url = 'https://www.murciasalud.es/archivo.php?id=470329'
download_file(tag, url)

# Pais Vasco - PV
tag = "PV"
remove_create(tag)
url = 'https://opendata.euskadi.eus/contenidos/ds_informes_estudios/covid_19_2020/opendata/situacion-epidemiologica.xlsx'
download_file(tag, url)
url = 'https://opendata.euskadi.eus/contenidos/ds_informes_estudios/covid_19_2020/opendata/datos-asistenciales.zip'
download_file(tag, url)
os.chdir("PV")
with zipfile.ZipFile("./datos-asistenciales.zip", 'r') as zip_ref:
  zip_ref.extractall("./")
os.chdir("..")

# Valencia - VC
tag = "VC"
remove_create(tag)
url = ''
download_file(tag, url)

# España
tag = "ES"
remove_create(tag)
url = 'https://cnecovid.isciii.es/covid19/resources/casos_diag_ccaadecl.csv'
download_file(tag, url)
url = 'https://cnecovid.isciii.es/covid19/resources/casos_hosp_uci_def_sexo_edad_provres.csv'
download_file(tag, url)

