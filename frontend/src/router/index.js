import Vue from "vue";
import Router from "vue-router";
import store from "@/store/index";

Vue.use(Router);

export const CHANGE_PASSWORD_PATH = "/account/password";

const router = new Router({
  mode: "history",
  base: '/',
  routes: [
    {
      path: "/",
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
      path: CHANGE_PASSWORD_PATH,
      component: () => import("@/views/ChangePasswordPage.vue"),
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
          path: "ResultFiles",
          component: () => import("@/components/modals/ResultFiles.vue"),
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
        {
          path: "workflow",
          component: () => import("@/components/adminComponents/WorkflowAdmin.vue")
        }
      ],
    },
    {
      path: "/tutorial",
      component: () => import("@/views/TutorialPage.vue"),
    },
  ],
});

// 초기·임시 비밀번호 상태에서는 비밀번호 변경 화면 외 이동 차단 (서버도 403으로 차단)
export function passwordChangeGuard(to, from, next) {
  const allowed = [CHANGE_PASSWORD_PATH, "/login"];
  if (store.getters.isLogin && store.getters.mustChangePassword && !allowed.includes(to.path)) {
    next(CHANGE_PASSWORD_PATH);
    return;
  }
  next();
}

router.beforeEach(passwordChangeGuard);

export default router;
