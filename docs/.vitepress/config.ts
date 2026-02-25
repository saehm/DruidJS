import { defineConfig } from "vitepress";
import typedocSidebar from "../api/typedoc-sidebar.json";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  markdown: {
    math: true,
  },
  title: "DruidJS",
  description: "Documentation",
  themeConfig: {
    siteTitle: "DruidJS",
    logo: "https://raw.githubusercontent.com/saehm/DruidJS/refs/heads/master/icon.svg",
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Showcases", link: "/showcase" },
      { text: "API Reference", link: "/api" },
    ],

    sidebar: [
      {
        text: "Showcases",
        link: "/showcase",
        collapsed: true,
        items: [
          { text: "Projections", link: "/showcase/projections" },
          { text: "Clustering", link: "/showcase/clustering" },
          { text: "Optimization", link: "/showcase/optimization" },
          { text: "Metrics", link: "/showcase/metrics" },
          { text: "Topology", link: "/showcase/topology" },
          { text: "Earth Mover", link: "/showcase/distribution" },
          { text: "Geospatial", link: "/showcase/geospatial" },
          { text: "Image Search", link: "/showcase/search" },
          { text: "Dendrograms", link: "/showcase/dendrogram" },
          { text: "Discovery", link: "/showcase/discovery" },
          { text: "K-Nearest Neighbors", link: "/showcase/knn" },
        ],
      },
      {
        text: "Dimensionality Reduction",
        collapsed: true,
        link: "/dimred",
        items: [
          { text: "PCA", link: "/dimred/pca" },
          {
            text: "FastMap",
            link: "/dimred/fastmap",
          },
          {
            text: "ISOMAP",
            link: "/dimred/isomap",
          },
          {
            text: "LDA",
            link: "/dimred/lda",
          },
          {
            text: "LLE",
            link: "/dimred/lle",
          },
          {
            text: "LSP",
            link: "/dimred/lsp",
          },
          {
            text: "LTSA",
            link: "/dimred/ltsa",
          },
          {
            text: "MDS",
            link: "/dimred/mds",
          },
          {
            text: "Sammon",
            link: "/dimred/sammon",
          },
          {
            text: "SMACOF",
            link: "/dimred/smacof",
          },
          {
            text: "SQDMDS",
            link: "/dimred/sqdmds",
          },
          {
            text: "t-SNE",
            link: "/dimred/tsne",
          },
          {
            text: "TopoMap",
            link: "/dimred/topomap",
          },
          {
            text: "TriMap",
            link: "/dimred/trimap",
          },
          {
            text: "UMAP",
            link: "/dimred/umap",
          },
        ],
      },
      {
        text: "API Reference",
        items: [...typedocSidebar],
      },
      { text: "LLMS", link: "/api/llms.txt" },
    ],

    socialLinks: [{ icon: "github", link: "https://github.com/saehm/DruidJS" }],

    search: {
      provider: "local",
    },
  },
});
