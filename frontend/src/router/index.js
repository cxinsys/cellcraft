import Vue from "vue";
import Router from "vue-router";

Vue.use(Router);

export default new Router({
  mode: "history",
  routes: [
    {
      path: "/",
      redirect: "/main",
    },
    {
      path: "/main",
      component: () => import("@/views/MainPage.vue"),
    },
    {
      path: "/login",
      component: () => import("@/views/LoginPage.vue"),
    },
    {
      path: "/signup",
      component: () => import("@/views/SignupPage.vue"),
    },
    {
      path: "/workflow",
      component: () => import("@/views/WorkFlowPage.vue"),
      children: [
        {
          path: "dataTable",
          component: () => import("@/components/modals/datatable.vue"),
        },
        {
          path: "file",
          component: () => import("@/components/modals/fileupload.vue"),
        },
        {
          path: "scatterPlot",
          component: () => import("@/components/modals/scatterPlot.vue"),
        },
        {
          path: "heatmap",
          component: () => import("@/components/modals/heatMap.vue"),
        },
        {
          path: "algorithm",
          component: () => import("@/components/modals/algorithm.vue"),
        },
      ],
    },
    {
      path: "/files",
      component: () => import("@/views/FilesPage.vue"),
    },
    {
      path: "/projects",
      component: () => import("@/views/ProjectPage.vue"),
    },
    {
      path: "/admin",
      component: () => import("@/views/AdminPage.vue"),
      children: [
        {
          path: "/",
          redirect: "/admin/user",
        },
        {
          path: "user",
          component: () => import("@/components/UserAdmin.vue"),
        },
        {
          path: "job",
          component: () => import("@/components/JobAdmin.vue"),
        },
        {
          path: "dataset",
          component: () => import("@/components/DatasetBoardAdmin.vue"),
        },
        {
          path: "algorithm",
          component: () => import("@/components/AlgorithmAdmin.vue"),
        },
      ],
    },
  ],
});
