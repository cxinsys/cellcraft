import Vue from 'vue'
import Router from 'vue-router'
// import CellcraftSidebar from '@/components/CellcraftSidebar'
// import CellcraftSidebar from './components/CellcraftSidebar'
import CellcraftSidebar from '../components/CellcraftSidebar'
import CellcraftHeader from '../components/CellcraftHeader'
Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/',
      name: 'CellcraftSidebar',
      component: CellcraftSidebar
    },
    {
      path: '/',
      name: 'CellcraftHeader',
      component: CellcraftHeader
    }
  ]
})
