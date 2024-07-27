import tkinter as tk
from tkinter import ttk, messagebox, Canvas, Scrollbar, Frame, Label, Entry, PanedWindow
import re
import matplotlib.pyplot as plt
import numpy as np
import random

# Função para preencher a região interna do polígono
# Funcionamento:
# A função `preenche_regiao_interna` é responsável por preencher a região interna de um polígono em uma matriz com a cor especificada.
# Ela percorre a matriz linha por linha, identifica os pontos que estão dentro do polígono com base nas bordas do polígono e 
# preenche esses pontos com a cor fornecida. A função utiliza a técnica de varredura horizontal para determinar quais pontos 
# estão dentro do polígono e garante que a área interna seja corretamente colorida.
def preenche_regiao_interna(matrix, cor):
    horizontal = []  # Lista para armazenar os pontos internos do polígono
    for i in range(matrix.shape[0]):  # Itera sobre as linhas da matriz
        inside_polygon = False  # Indicador se o ponto está dentro do polígono
        horizontal_linha = []  # Lista temporária para armazenar os pontos da linha atual
        counter = 0  # Contador de bordas do polígono encontradas
        for j in range(matrix.shape[1]):  # Itera sobre as colunas da matriz
            # Verifica se o ponto atual faz parte da borda do polígono
            if all(matrix[i][j] == cor):
                # Verifica se não é o último ponto da linha
                if j != matrix.shape[1] - 1:
                    # Verifica se o próximo ponto não faz parte da borda do polígono
                    if all(matrix[i][j + 1] != cor):
                        # Alterna o estado do indicador (dentro/fora do polígono)
                        inside_polygon = not inside_polygon
                        counter += 1  # Incrementa o contador de bordas
                else:
                    inside_polygon = not inside_polygon
                    counter += 1
            elif inside_polygon:
                horizontal_linha.append((i, j))  # Adiciona o ponto à lista temporária
        if counter > 1:
            # Adiciona os pontos internos da linha atual à lista geral
            for elemento in horizontal_linha:
                horizontal.append(elemento)
    # Preenche os pontos internos do polígono com a cor especificada
    for point in horizontal:
        matrix[point[0]][point[1]] = cor

# Função para normalizar as coordenadas
# Funcionamento:
# A função `nova_coord` é responsável por normalizar as coordenadas de um intervalo de [-1, 1]
# para o intervalo correspondente dentro de uma matriz com largura W e altura H.
# Essa normalização é necessária para ajustar coordenadas de um sistema de referência para o sistema
# de coordenadas da matriz, garantindo que as coordenadas se ajustem corretamente ao tamanho da imagem.
def nova_coord(old_x, old_y, W, H):
    # Normaliza a coordenada x
    # Converte old_x do intervalo [-1, 1] para o intervalo [0, W]
    new_x = (old_x - (-1)) * W / 2
    
    # Normaliza a coordenada y
    # Converte old_y do intervalo [-1, 1] para o intervalo [0, H]
    new_y = (old_y - (-1)) * H / 2
    
    # Retorna as novas coordenadas normalizadas
    return [new_x, new_y]

# Função para produzir um fragmento na matriz
# Funcionamento:
# A função `produz_fragmento` é responsável por definir a cor de um ponto específico na matriz.
def produz_fragmento(x, y, matrix, cor):
    Xm = round(x)  # Arredonda a coordenada x para o valor inteiro mais próximo
    Ym = round(y)  # Arredonda a coordenada y para o valor inteiro mais próximo
    
    # Garante que as coordenadas estejam dentro dos limites da matriz
    Xm = min(max(0, Xm), len(matrix[0]) - 1)  # Ajusta Xm para que esteja entre 0 e largura da matriz - 1
    Ym = min(max(0, Ym), len(matrix) - 1)      # Ajusta Ym para que esteja entre 0 e altura da matriz - 1
    
    try:
        matrix[Ym][Xm] = cor  # Aplica a cor no ponto especificado da matriz
    except IndexError:
        return  # Em caso de erro de índice, simplesmente retorna sem fazer alterações

# Função de rasterização de uma reta
# Overview:
# A função `algoritmo_rasterizacao` implementa o algoritmo de rasterização de linhas para desenhar uma linha reta entre dois pontos em uma matriz.
# Utilizando o conceito de coordenadas normalizadas e o algoritmo de Bresenham, a função determina quais pontos devem ser coloridos para 
# representar a linha na matriz de pixels. Ela faz ajustes finos nos pontos para garantir que a linha seja desenhada corretamente no espaço 2D da matriz.
#
# Parâmetros:
# X1 (int): Coordenada x do ponto inicial da linha.
# Y1 (int): Coordenada y do ponto inicial da linha.
# X2 (int): Coordenada x do ponto final da linha.
# Y2 (int): Coordenada y do ponto final da linha.
# H (int): Altura da matriz de pixels onde a linha será desenhada.
# W (int): Largura da matriz de pixels onde a linha será desenhada.
# matrix (numpy.ndarray): Matriz 2D (ou 3D se colorida) que representa a imagem onde a linha será desenhada.
# cor (list): Lista com três valores (R, G, B) que representam a cor da linha a ser desenhada.

def algoritmo_rasterizacao(X1, Y1, X2, Y2, H, W, matrix, cor):
    # Normaliza as coordenadas dos pontos de início e fim da linha
    X1, Y1 = nova_coord(X1, Y1, H, W)
    X2, Y2 = nova_coord(X2, Y2, H, W)
    # Calcula as diferenças nas coordenadas x e y
    dx = X2 - X1
    dy = Y2 - Y1
    # Calcula o coeficiente angular da linha (slope)
    m = dy / dx if dx != 0 else 0
    # Inicializa as coordenadas do ponto atual com as coordenadas de início
    x, y = X1, Y1
    # Calcula o intercepto da linha
    b = y - m * x

    # Função auxiliar para desenhar o ponto baseado na coordenada x
    def func_plot_x(x, matrix):
        y = m * x + b  # Calcula a coordenada y correspondente
        produz_fragmento(x, y, matrix, cor)  # Desenha o ponto na matriz

    # Função auxiliar para desenhar o ponto baseado na coordenada y
    def func_plot_y(y, matrix):
        produz_fragmento(x, y, matrix, cor)  # Desenha o ponto na matriz

    # Define o incremento para o ajuste fino da linha
    incremento = 0.1
    
    # Desenha a linha da esquerda para a direita
    if x < X2:
        while x < X2:
            x += incremento
            func_plot_x(x, matrix)
    
    # Desenha a linha da direita para a esquerda
    elif x == X2:
        if y < Y2:
            while y < Y2:
                y += incremento
                func_plot_y(y, matrix)
        else:
            while y > Y2:
                y -= incremento
                func_plot_y(y, matrix)
    
    # Desenha a linha da direita para a esquerda
    else:
        while x > X2:
            x -= incremento
            func_plot_x(x, matrix)

# Função para checar se o range é válido
# Funcionamento:
# A função `checa_validade_do_range` verifica se as coordenadas fornecidas estão dentro dos limites permitidos de -1 a 1.
# Se alguma coordenada estiver fora desses limites, a função exibe uma mensagem de erro detalhando quais coordenadas estão inválidas,
# e para qual reta (ou polígono) elas pertencem. Esta função é útil para garantir que as coordenadas fornecidas para rasterização
# estejam dentro do intervalo esperado e evitem erros de processamento.
def checa_validade_do_range(input_list, reta_atual):
    tem_erro = False  # Flag para indicar se houve algum erro
    fora_dos_limites = []  # Lista para armazenar coordenadas que estão fora dos limites permitidos

    # Itera sobre cada coordenada na lista fornecida
    for index in range(len(input_list)):
        # Verifica se a coordenada está fora do intervalo permitido [-1, 1]
        if -1 > float(input_list[index]) or float(input_list[index]) > 1:
            tem_erro = True  # Define a flag de erro como True
            fora_dos_limites.append(input_list[index])  # Adiciona a coordenada fora dos limites à lista

    # Se houve algum erro, exibe uma mensagem de erro
    if tem_erro:
        messagebox.showerror(
            title='Erro',
            message=f"Coordenada(s) {fora_dos_limites} da reta {reta_atual} fora dos limites impostos (-1 a 1)!"
        )

# Função para desenhar arestas
# Funcionamento: A função `produz_arestas` desenha as arestas de um polígono em uma matriz de pixels usando o algoritmo de rasterização de linhas.
# Para cada aresta do polígono, a função chama o algoritmo de rasterização para desenhar a linha correspondente na matriz.
# Isso resulta na representação visual das arestas do polígono na matriz, que pode ser usada para visualização ou processamento posterior.

def produz_arestas(array_pontos, array_arestas, H, W, matrix, cor):
    print(f"array pontos {array_pontos} array arestas {array_arestas} cor {cor}")  # Imprime os pontos, arestas e cor para depuração
    # Itera sobre cada aresta no array de arestas
    for aresta in array_arestas:
        primeiro_ponto, segundo_ponto = aresta  # Obtém os índices dos pontos que formam a aresta
        origem_x, origem_y = array_pontos[primeiro_ponto]  # Obtém as coordenadas do ponto inicial da aresta
        destino_x, destino_y = array_pontos[segundo_ponto]  # Obtém as coordenadas do ponto final da aresta

        print(f"array_pontos[primeiro_ponto] {array_pontos[primeiro_ponto]} array_pontos[segundo_ponto] {array_pontos[segundo_ponto]}")  # Imprime as coordenadas dos pontos para depuração

        # Chama o algoritmo de rasterização para desenhar a linha (aresta) entre os dois pontos
        algoritmo_rasterizacao(origem_x, origem_y, destino_x, destino_y, H, W, matrix, cor)

# Função para desenhar polígonos
# Funcionamento: A função `desenha_poligono` desenha polígonos em uma matriz de pixels, preenchendo as arestas e a região interna do polígono.
# Para cada polígono no array de polígonos, a função gera as arestas e chama a função `produz_arestas` para desenhá-las.
# Após desenhar as arestas, a função `preenche_regiao_interna` é chamada para preencher a área interna do polígono com a cor especificada.

def desenha_poligono(array_pontos, array_poligonos, H, W, matrix, cor):
    # Itera sobre cada polígono no array de polígonos
    for poligono in array_poligonos:
        # Cria uma lista de arestas para o polígono atual
        array_arestas = [[i, (i + 1) % len(poligono)] for i in range(len(poligono))]
        # Desenha as arestas do polígono usando a função `produz_arestas`
        produz_arestas([array_pontos[i] for i in poligono], array_arestas, H, W, matrix, cor)
        # Preenche a região interna do polígono com a cor especificada
        preenche_regiao_interna(matrix, cor)

# Função para rasterizar retas
def rasterizacao_retas():


    # A função `produz_retas_continuas` desenha retas contínuas em um gráfico usando Matplotlib.
    # Ela itera sobre uma lista de pontos e suas respectivas cores para traçar linhas entre os pontos consecutivos.
    # As retas são desenhadas com base nos valores de cor especificados e são exibidas em um subplot da figura Matplotlib.
    def produz_retas_continuas(array_pontos, cores):
        # Itera sobre a lista de pontos e suas respectivas cores
        for index in range(len(array_pontos)):
            # Obtém a cor para a reta atual
            cor = cores[index]
            # Define uma lista de arestas conectando pontos consecutivos
            array_retas = [[0, 1]]
            # Itera sobre cada aresta na lista de arestas
            for aresta in array_retas:
                # Extrai os índices dos pontos que formam a aresta
                primeiro_ponto, segundo_ponto = aresta
                # Obtém as coordenadas dos pontos de origem e destino da reta
                origem_x, origem_y = array_pontos[index][primeiro_ponto]
                destino_x, destino_y = array_pontos[index][segundo_ponto]
                # Adiciona a reta ao gráfico atual usando Matplotlib
                plt.subplot(2, 3, 6)  # Seleciona o subplot onde a reta será desenhada
                plt.xlim([-1, 1])     # Define o intervalo do eixo x
                plt.ylim([-1, 1])     # Define o intervalo do eixo y
                plt.plot([origem_x, destino_x], [origem_y, destino_y],
                        color=(cor[0] / 255, cor[1] / 255, cor[2] / 255))  # Plota a reta com a cor especificada


    # A função `desenha_reta` é responsável por desenhar retas em diferentes resoluções e exibir essas imagens usando Matplotlib.
    # Ela cria uma matriz para cada resolução especificada, desenha as retas nas matrizes e, em seguida, exibe essas matrizes como imagens.
    # Além disso, também plota as retas contínuas em um gráfico para comparação visual.
    def desenha_reta(array_pontos):
        # Define uma lista de resoluções para as imagens a serem geradas
        resolucoes = [[100, 100], [300, 300], [600, 600], [800, 600], [1920, 1080]]
        # Gera uma lista de cores aleatórias para cada ponto
        cores = [[random.randint(0, 255) for _ in range(3)] for _ in range(len(array_pontos))]
        
        # Itera sobre cada resolução especificada
        for resolucao in resolucoes:
            # Cria uma matriz de zeros (imagem preta) com a resolução especificada
            matrix = np.zeros((resolucao[1], resolucao[0], 3), dtype=np.uint8)
            
            # Itera sobre cada conjunto de pontos na lista de pontos
            for index in range(len(array_pontos)):
                # Define a lista de arestas a serem desenhadas (aqui conectando pontos consecutivos)
                array_retas = [[0, 1]]
                # Chama a função para desenhar as arestas na matriz
                produz_arestas(array_pontos[index], array_retas, resolucao[0], resolucao[1], matrix, cores[index])
            
            # Seleciona o subplot na posição apropriada para exibir a imagem
            plt.subplot(2, 3, resolucoes.index(resolucao) + 1)
            # Exibe a matriz como uma imagem
            plt.imshow(matrix.astype("uint8"))
            # Inverte o eixo y para a orientação correta da imagem
            plt.gca().invert_yaxis()
        
        # Plota as retas contínuas em um gráfico separado para comparação
        produz_retas_continuas(array_pontos, cores)
        # Exibe todos os gráficos e imagens gerados
        plt.show()

    # A função `pegar_retas` é responsável por coletar as coordenadas das retas a partir da interface gráfica, validar as coordenadas, e então chamar a função para desenhar essas retas.
    # Ela lida com a extração e validação de coordenadas de entrada, verifica se as coordenadas estão dentro dos limites permitidos, e executa o desenho das retas se os dados forem válidos.

    def pegar_retas():
        coords = []  # Lista para armazenar as coordenadas das retas
        tem_erro = False  # Flag para indicar se houve algum erro na coleta de dados
        count = 0  # Contador para rastrear o número da reta atual
        
        # Itera sobre cada entrada de reta na lista lbllst
        for reta in lbllst:
            count += 1  # Incrementa o contador de retas
            # Extrai números das coordenadas fornecidas usando expressão regular
            input_list = re.findall(r'-?\d+\.?\d*', reta[1].get())
            
            # Verifica se a lista de coordenadas está vazia
            if len(input_list) == 0:
                tem_erro = True  # Define a flag de erro como True
                messagebox.showerror(title='Erro', message=f"Não foi possível identificar coordenadas na reta {count}!")  # Exibe uma mensagem de erro
            # Verifica se a lista de coordenadas tem menos de 4 elementos
            elif len(input_list) < 4:
                tem_erro = True  # Define a flag de erro como True
                messagebox.showerror(title='Erro', message=f"Na reta {count} você tem pontos insuficientes para printar a reta!")  # Exibe uma mensagem de erro
            # Verifica se a lista de coordenadas tem mais de 4 elementos
            elif len(input_list) > 4:
                tem_erro = True  # Define a flag de erro como True
                messagebox.showerror(title='Erro', message=f"Na reta {count} você tem pontos a mais do que os necessários para printar a reta!")  # Exibe uma mensagem de erro
            
            # Se não houve erros na coleta de dados
            if not tem_erro:
                # Verifica se as coordenadas estão dentro dos limites permitidos
                checa_validade_do_range(input_list, count)
            
            # Se não houve erros na coleta de dados e as coordenadas são válidas
            if not tem_erro:
                # Converte as coordenadas para floats e as adiciona à lista de coordenadas
                coords.append([[float(input_list[0]), float(input_list[1])], [float(input_list[2]), float(input_list[3])]])
        
        # Se não houve erros durante a coleta e validação de dados
        if not tem_erro:
            # Chama a função para desenhar as retas com as coordenadas coletadas
            desenha_reta(coords)

    # Interface de usuário para retas
    def abrir_retas():
        reta_window = tk.Toplevel()
        reta_window.title("Rasterização de Retas")
        reta_window.geometry("500x400+750+300")
        reta_window.resizable(False, False)

        # Valores de teste
        valores_padrao_retas = [
            "0 -0.5 0.5 0.5",
            "0 0.3 0.8 0.8"
        ]



        panel1 = PanedWindow(reta_window, bd=4, relief='raised', orient='horizontal')
        panel1.pack(fill='both', expand=1)

        left_label = tk.Label(panel1)
        panel1.add(left_label)

        panel2 = PanedWindow(panel1, orient='vertical', bd=4, relief='raised')
        panel1.add(panel2)

        right_label = tk.Label(panel2, text="INPUTS")
        panel2.add(right_label)

        scrollable_frame = tk.Frame(right_label)
        scrollable_frame.pack(fill='both', expand=1)

        scrollbar = ttk.Scrollbar(scrollable_frame, orient='vertical')
        scrollbar.pack(side='right', fill='y')

        canvas = Canvas(scrollable_frame, scrollregion=(0, 0, 0, 8000), highlightthickness=0)
        canvas.pack(side='left', fill='both', expand=1)
        scrollbar.config(command=canvas.yview)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        lbllst.clear()

        for i in range(len(valores_padrao_retas)):
            frame = tk.Frame(content_frame)
            frame.pack(pady=10)

            label = tk.Label(frame, text=f"Reta {i + 1} (x1 y1 x2 y2):")
            label.pack(side='left')

            retas = tk.Entry(frame)
            retas.pack(side='left', padx=10)
            retas.insert(0, valores_padrao_retas[i])  # Adiciona valor padrão

            lbllst.append([label, retas])

        btn_sair = tk.Button(reta_window, text="Sair", command=reta_window.destroy)
        btn_sair.pack(pady=10)

        btn_abrir = tk.Button(reta_window, text="Rasterizar", command=pegar_retas)
        btn_abrir.pack(pady=10)

        reta_window.mainloop()

    lbllst = []
    abrir_retas()

# Função para rasterizar polígonos
def rasterizacao_poligonos():
    # Função para desenhar polígonos em uma matriz
    # Overview:
    # A função `desenha_poligonos` é responsável por desenhar uma série de polígonos em uma matriz de pixels e exibir os resultados usando matplotlib.
    # Para cada polígono, a função gera uma matriz de pixels coloridos onde o polígono é desenhado. Ela usa a função `desenha_poligono` para realizar a
    # rasterização de cada polígono e cria imagens em diferentes resoluções, exibindo-as em um grid de subplots.
    #
    # Parâmetros:
    # array_pontos (list of lists): Lista contendo as coordenadas dos vértices dos polígonos. Cada sublista representa os pontos (vértices) de um polígono.

    def desenha_poligonos(array_pontos):
        # Define as resoluções para a matriz de pixels
        resolucoes = [[100, 100], [300, 300], [600, 600], [800, 600], [1920, 1080]]
        # Gera cores aleatórias para cada polígono
        cores = [[random.randint(0, 255) for _ in range(3)] for _ in range(len(array_pontos))]
        
        # Itera sobre cada resolução definida
        for resolucao in resolucoes:
            # Cria uma matriz de pixels com a resolução especificada
            matrix = np.zeros((resolucao[1], resolucao[0], 3), dtype=np.uint8)
            
            # Itera sobre cada polígono definido em array_pontos
            for index in range(len(array_pontos)):
                # Cria um array de polígonos com um único polígono para a iteração atual
                array_poligonos = [list(range(len(array_pontos[index])))]
                # Desenha o polígono na matriz usando a função desenha_poligono
                desenha_poligono(array_pontos[index], array_poligonos, resolucao[0], resolucao[1], matrix, cores[index])
            
            # Exibe a matriz de pixels como uma imagem usando matplotlib
            plt.subplot(2, 3, resolucoes.index(resolucao) + 1)
            plt.imshow(matrix.astype("uint8"))
            plt.gca().invert_yaxis()
        
        # Mostra todas as imagens em um grid de subplots
        plt.show()

    # Pega as coordenadas dos polígonos digitados pelo usuário
    def pegar_poligonos():
        coords = []
        tem_erro = False
        count = 0
        for poligono in lbllst:
            count += 1
            input_list = re.findall(r'-?\d+\.?\d*', poligono[1].get())
            print("lista", input_list)
            if len(input_list) == 0:
                tem_erro = True
                messagebox.showerror(title='Erro', message=f"Não foi possível identificar coordenadas do polígono {count}!")
            else:
                checa_validade_do_range(input_list, count)
                array_pontos = []
                for index in range(0, len(input_list), 2):
                    array_pontos.append([float(input_list[index]), float(input_list[index + 1])])
                coords.append(array_pontos)
        if not tem_erro:
            desenha_poligonos(coords)

    # Interface de usuário para polígonos
    def abrir_poligonos():
        poligono_window = tk.Toplevel()
        poligono_window.title("Rasterização de Polígonos")
        poligono_window.geometry("500x400+750+300")
        poligono_window.resizable(False, False)

        # Valores padrão para polígonos
        valores_padrao_poligonos = [
            "0 0 0.5 0 0.25 0.5",  # Triângulo
            "0 -0.5 0.5 -0.5 0.5 0 0 0",  # Quadrado
            "0 0 0.5 0.5 -0.5 0.5"  # Losango
        ]


        panel1 = PanedWindow(poligono_window, bd=4, relief='raised', orient='horizontal')
        panel1.pack(fill='both', expand=1)

        left_label = tk.Label(panel1)
        panel1.add(left_label)

        panel2 = PanedWindow(panel1, orient='vertical', bd=4, relief='raised')
        panel1.add(panel2)

        right_label = tk.Label(panel2, text="INPUTS")
        panel2.add(right_label)

        scrollable_frame = tk.Frame(right_label)
        scrollable_frame.pack(fill='both', expand=1)

        scrollbar = ttk.Scrollbar(scrollable_frame, orient='vertical')
        scrollbar.pack(side='right', fill='y')

        canvas = Canvas(scrollable_frame, scrollregion=(0, 0, 0, 8000), highlightthickness=0)
        canvas.pack(side='left', fill='both', expand=1)
        scrollbar.config(command=canvas.yview)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        lbllst.clear()

        for i in range(len(valores_padrao_poligonos)):
            frame = tk.Frame(content_frame)
            frame.pack(pady=10)

            label = tk.Label(frame, text=f"Polígono {i + 1} (x1 y1 x2 y2 ...):")
            label.pack(side='left')

            poligono = tk.Entry(frame)
            poligono.pack(side='left', padx=10)
            poligono.insert(0, valores_padrao_poligonos[i])  # Adiciona valor padrão

            lbllst.append([label, poligono])

        btn_sair = tk.Button(poligono_window, text="Sair", command=poligono_window.destroy)
        btn_sair.pack(pady=10)

        btn_abrir = tk.Button(poligono_window, text="Rasterizar", command=pegar_poligonos)
        btn_abrir.pack(pady=10)

        poligono_window.mainloop()

    lbllst = []
    abrir_poligonos()

# Função para rasterizar curvas
def rasterizacao_curvas():
    # Função para desenhar curvas usando Hermite
    def desenha_curvas(array_curvas):
        # Define as resoluções para a matriz de pixels
        resolucoes = [[100, 100], [300, 300], [600, 600], [800, 600], [1920, 1080]]
        
        # Itera sobre cada resolução definida
        for resolucao in resolucoes:
            # Cria uma matriz de pixels com a resolução especificada
            matrix = np.zeros((resolucao[1], resolucao[0], 3), dtype=np.uint8)
            
            # Itera sobre cada curva definida em array_curvas
            for curva in array_curvas:
                P1, P2, T1, T2, cor = curva
                # Desenha a curva na matriz usando a função de rasterização específica
                desenha_curva_hermite(P1, P2, T1, T2, resolucao[0], resolucao[1], matrix, cor)
            
            # Exibe a matriz de pixels como uma imagem usando matplotlib
            plt.subplot(2, 3, resolucoes.index(resolucao) + 1)
            plt.imshow(matrix.astype("uint8"))
            plt.gca().invert_yaxis()
        
        # Mostra todas as imagens em um grid de subplots
        plt.show()

    # Função para desenhar uma curva Hermite na matriz
    def desenha_curva_hermite(P1, P2, T1, T2, W, H, matrix, cor):
        # Gera 100 pontos igualmente espaçados entre 0 e 1
        t = np.linspace(0, 1, num=100)
        
        # Calcula as funções base de Hermite
        H0 = (2 * t**3 - 3 * t**2 + 1)  # Função base H0
        H1 = (-2 * t**3 + 3 * t**2)     # Função base H1
        H2 = (t**3 - 2 * t**2 + t)      # Função base H2
        H3 = (t**3 - t**2)              # Função base H3
        
        # Calcula as coordenadas x e y da curva usando a fórmula de Hermite
        x = H0 * P1[0] + H1 * P2[0] + H2 * T1[0] + H3 * T2[0]
        y = H0 * P1[1] + H1 * P2[1] + H2 * T1[1] + H3 * T2[1]
        
        # Desenha a curva ponto a ponto
        for i in range(len(x)):
            # Converte as coordenadas (x[i], y[i]) para o sistema de coordenadas da matriz
            x_norm, y_norm = nova_coord(x[i], y[i], W, H)
            
            # Desenha o ponto na matriz na posição (x_norm, y_norm) com a cor especificada
            produz_fragmento(x_norm, y_norm, matrix, cor)

    # Função para coletar as curvas da interface gráfica e validar as entradas
    def pegar_curvas():
        coords = []
        tem_erro = False
        count = 0
        
        for curva in lbllst:
            count += 1
            input_list = re.findall(r'-?\d+\.?\d*', curva[1].get())
            
            if len(input_list) != 8:
                tem_erro = True
                messagebox.showerror(title='Erro', message=f"Na curva {count}, você deve fornecer 8 valores (P1x P1y P2x P2y T1x T1y T2x T2y)!")
            else:
                checa_validade_do_range(input_list, count)
                P1 = [float(input_list[0]), float(input_list[1])]
                P2 = [float(input_list[2]), float(input_list[3])]
                T1 = [float(input_list[4]), float(input_list[5])]
                T2 = [float(input_list[6]), float(input_list[7])]
                cor = [random.randint(0, 255) for _ in range(3)]
                coords.append([P1, P2, T1, T2, cor])
        
        if not tem_erro:
            desenha_curvas(coords)

    # Interface de usuário para curvas
    def abrir_curvas():
        curva_window = tk.Toplevel()
        curva_window.title("Rasterização de Curvas")
        curva_window.geometry("500x400+750+300")
        curva_window.resizable(False, False)

        # Valores padrão para curvas
        valores_padrao_curvas = [
            "0 0 0.5 0 0 0.5 0.5 0",
            "0 0 0.3 0.3 0.1 0.5 0.4 0.5"
        ]

        panel1 = PanedWindow(curva_window, bd=4, relief='raised', orient='horizontal')
        panel1.pack(fill='both', expand=1)

        left_label = tk.Label(panel1)
        panel1.add(left_label)

        panel2 = PanedWindow(panel1, orient='vertical', bd=4, relief='raised')
        panel1.add(panel2)

        right_label = tk.Label(panel2, text="INPUTS")
        panel2.add(right_label)

        scrollable_frame = tk.Frame(right_label)
        scrollable_frame.pack(fill='both', expand=1)

        scrollbar = ttk.Scrollbar(scrollable_frame, orient='vertical')
        scrollbar.pack(side='right', fill='y')

        canvas = Canvas(scrollable_frame, scrollregion=(0, 0, 0, 8000), highlightthickness=0)
        canvas.pack(side='left', fill='both', expand=1)
        scrollbar.config(command=canvas.yview)

        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor='nw')

        lbllst.clear()

        for i in range(len(valores_padrao_curvas)):
            frame = tk.Frame(content_frame)
            frame.pack(pady=10)

            label = tk.Label(frame, text=f"Curva {i + 1} (P1x P1y P2x P2y T1x T1y T2x T2y):")
            label.pack(side='left')

            curvas = tk.Entry(frame)
            curvas.pack(side='left', padx=10)
            curvas.insert(0, valores_padrao_curvas[i])  # Adiciona valor padrão

            lbllst.append([label, curvas])

        btn_sair = tk.Button(curva_window, text="Sair", command=curva_window.destroy)
        btn_sair.pack(pady=10)

        btn_abrir = tk.Button(curva_window, text="Rasterizar", command=pegar_curvas)
        btn_abrir.pack(pady=10)

        curva_window.mainloop()

    lbllst = []
    abrir_curvas()

# Função principal para iniciar a aplicação
def main():
    root = tk.Tk()
    root.title("Rasterização de Polígonos e Retas")
    root.geometry("300x300")

    btn_rasterizar_retas = tk.Button(root, text="Rasterizar Retas", command=rasterizacao_retas)
    btn_rasterizar_retas.pack(pady=10)

    btn_rasterizar_poligonos = tk.Button(root, text="Rasterizar Polígonos", command=rasterizacao_poligonos)
    btn_rasterizar_poligonos.pack(pady=10)

    btn_rasterizar_curvas = tk.Button(root, text="Rasterizar Curvas", command=rasterizacao_curvas)
    btn_rasterizar_curvas.pack(pady=10)

    root.mainloop()

if __name__ == "__main__":
    main()