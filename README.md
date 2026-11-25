# Análise Espaço-Temporal do Índice Espectral CO₂Flux no Perímetro Urbano de Barreiras-BA (2016–2025)

Repositório contendo o código-fonte em JavaScript desenvolvido para a plataforma **Google Earth Engine (GEE)**, referente ao Trabalho de Conclusão de Curso (TCC) apresentado ao Curso de Bacharelado em Ciências Biológicas da **Universidade Federal do Oeste da Bahia (UFOB)**.

- **Autora:** Leci Silva França
- **Orientador:** Prof. Dr. Elvis Bergue Mariz Moreira
- **Ano:** 2026

## Sobre o Projeto

O script realiza o processamento digital de imagens orbitais em ambiente de computação em nuvem para analisar a dinâmica espaço-temporal do potencial de assimilação de carbono da vegetação intraurbana de Barreiras, Bahia, ao longo do decênio de **2016 a 2025**, totalizando **1.065 cenas orbitais válidas** da constelação Sentinel-2.

### Principais etapas executadas no script:
1. **Aquisição e Fusão de Dados (*Data Fusion*):** Carregamento da coleção harmonizada Sentinel-2 (`COPERNICUS/S2_HARMONIZED`, resolução espacial de 10 m) e vinculação à coleção de probabilidade de nuvens (`COPERNICUS/S2_CLOUD_PROBABILITY`) via `ee.Join.saveFirst()`.
2. **Filtragem de Nuvens e Efeitos de Borda:** Mascaramento probabilístico (< 40%) e remoção de pixels ruidosos nas bordas das órbitas utilizando as bandas B8A e B9.
3. **Modelagem Espectral (`CO2_Index`):** Cálculo anual por mediana dos índices **NDVI**, **PRI adaptado** (`PRI_proxy`), **PRI escalonado** (`sPRI`) e do modelo integrativo **CO₂Flux** (`NDVI × sPRI`).
4. **Análise Estatística Não Paramétrica:** Execução pixel a pixel do teste de **Mann-Kendall** (direção e significância estatística com α = 0,05) e do **Estimador de Declive de Sen** (*Sen's Slope* - magnitude anual).
5. **Quantificação Territorial e Exportação:** Cálculo das áreas em hectares (ha) e porcentagem (%) para as cinco classes de tendência e exportação dos produtos tabulares (`.csv`) e matriciais (`.tif`).

## Como executar
1. Acesse o Google Earth Engine Code Editor (https://code.earthengine.google.com/).
2. Copie o conteúdo do arquivo `script_co2flux_barreiras.js` deste repositório e cole no editor.
3. Certifique-se de ter acesso ao vetor de delimitação da área de estudo (`projects/elvisbergmm/assets/Limites_Bairros`) ou substitua pelo seu próprio *Asset* geométrico.
4. Clique em **Run** para executar o processamento e acompanhar os gráficos, a contagem de cenas e as tarefas de exportação nas abas *Console* e *Tasks*.
## Como citar

Caso utilize este código ou adapte a metodologia em suas pesquisas, por favor, cite o repositório original:

> FRANÇA, Leci Silva. **CO2Flux-Barreiras-BA-GEE**: análise espaço-temporal do índice espectral CO2Flux no perímetro urbano de Barreiras-BA (2016-2025). Barreiras: GitHub, 2026. Repositório de código-fonte em JavaScript (Google Earth Engine). Disponível em: https://github.com/LeciSilva/CO2Flux-Barreiras-BA-GEE.
