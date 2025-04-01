import Vue from "vue";
import Router from "vue-router";

Vue.use(Router);

export default new Router({
  mode: "history",
  base: '/cellcraft',
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
          path: "InputFile",
          component: () => import("@/components/modals/InputFile.vue"),
        },
        {
          path: "Datatable",
          component: () => import("@/components/modals/DataTable.vue"),
        },
        {
          path: "ScatterPlot",
          component: () => import("@/components/modals/ScatterPlot.vue"),
        },
        {
          path: "Algorithm",
          component: () => import("@/components/modals/Algorithm.vue"),
        },
        {
          path: "ResultFile",
          component: () => import("@/components/modals/ResultFile.vue"),
        },
        {
          path: "Visualization",
          component: () => import("@/components/modals/Visualization.vue"),
        }
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
      path: "/datasets",
      component: () => import("@/views/DatasetsPage.vue"),
    },
    {
      path: "/plugins",
      component: () => import("@/views/PluginsPage.vue"),
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
          component: () => import("@/components/adminComponents/UserAdmin.vue"),
        },
        {
          path: "job",
          component: () => import("@/components/adminComponents/TaskAdmin.vue"),
        },
        {
          path: "dataset",
          component: () => import("@/components/adminComponents/FileAdmin.vue"),
        },
        {
          path: "algorithm",
          component: () => import("@/components/adminComponents/PluginAdmin.vue"),
        },
      ],
    },
    {
      path: "/tutorial",
      component: () => import("@/views/TutorialPage.vue"),
    },
  ],
});
