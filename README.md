# NetAIR
R Network Analysis Interface based in plumber API

## Instalação
- Instale o R.
- Instale o RStudio.
- Na página deste repositório, clique no botão verde *Clone or Download*, em seguida em *Download ZIP*.
- Extraia os arquivos em uma pasta de sua preferência.
- Abra a pasta principal e clique em NetAIR (abrirá uma nova sessão no RStudio).
- Abra o arquivo **plumber.R**, instale os *packages* indicados, execute o código (botão Run API, abrirá uma janela no seu navegador).

## Análise
- Substitua */__swagger__/* na barra de endereço por */start* e pressione **enter**. Exemplo: 999.999.99/start.
- Coloque seu banco de dados em formato *.csv* na pasta *data* (contida na pasta principal que você baixou deste repositório).
  - Seu banco de dados não deve conter *missing values*.
  - Na segunda linha você deve especificar o tipo de variável (g = Gaussiana, p = Poisson, c = Categórica).
- Digite o nome da variável dicotômica pela qual você quer dividir seu banco de dados e gerar as redes. 
- Digite o número de permutações.
- Clique em **Submeter** e espere alguns minutos (por exemplo, 1000 permutações pode levar horas).
- Quando a análise estiver concluída, clique em *gerar relatório*.
