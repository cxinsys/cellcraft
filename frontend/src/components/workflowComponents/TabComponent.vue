<template>
    <VueDragResize contentClass="content-component" v-if="tabList.length != 0 && isTabView" :isActive="true" :x="600"
        :y="64" :w="880" :h="672" :minw="820" :minh="540" :stickSize="14" :sticks="['tl', 'ml', 'tr', 'bl', 'br']">
        <ul class="content-tab" v-if="tabList.length != 0 && isTabView">
            <li class="tab__item" v-for="(tab, idx) in tabList" :key="idx"
                v-bind:class="{ currentTab: currentTab === idx }" @click="tabClick(idx)">
                <div class="tab__name">
                    <img class="tab__icon" :src="tab.img" />
                    <p class="tab__text">
                        {{ tab.title || tab.name }}
                    </p>
                </div>
                <img class="tab__close" @click="removeTab(tab.id)" src="@/assets/close.png" />
            </li>
        </ul>
        <div class="tab__hide" @click="$emit('update:isTabView', false)"></div>
        <div class="content-view" v-if="tabList.length != 0 && isTabView" @mousedown.stop>
            <router-view :key="$route.fullPath"></router-view>
        </div>
    </VueDragResize>
</template>

<script>
import VueDragResize from "vue-drag-resize";
export default {
    name: 'TabComponent',
    props: {
        initialTabList: {
            type: Array,
            default: () => []
        },
        isTabView: {
            type: Boolean,
            required: true
        },
        currentWorkflowId: {
            type: Number,
            required: true
        }
    },
    components: {
        VueDragResize,
    },
    data() {
        return {
            tabList: [...this.initialTabList],
            currentTab: 0
        }
    },
    methods: {
        createTab(node) {
            const newTab = {
                id: node.id,
                name: node.name,
                title: node.data.title || node.name,
                img: require(`@/assets/${node.name}.png`)
            };
            this.tabList.push(newTab);
            this.currentTab = this.tabList.length - 1;
            this.componentChange(newTab);
        },
        removeTab(id) {
            const index = this.tabList.findIndex(tab => tab.id === id);
            if (index !== -1) {
                this.tabList.splice(index, 1);
                if (this.currentTab >= index) {
                    // this.currentTab = this.tabList.length - 1;
                    this.currentTab = this.currentTab === 0 ? 0 : this.currentTab - 1;
                }
                if (this.tabList.length > 0) {
                    this.componentChange(this.tabList[this.currentTab]);
                }
            }
        },
        adjustCurrentTab(node) {
            const index = this.tabList.findIndex(tab => tab.id === node.id);
            if (index !== -1) {
                this.currentTab = index;
            } else {
                this.createTab(node);
            }
            if (this.tabList.length > 0) {
                this.componentChange(this.tabList[this.currentTab]);
            }
        },
        tabClick(idx) {
            this.currentTab = idx;
            this.componentChange(this.tabList[idx]);
        },
        componentChange(tab) {
            let newPath = `/workflow/${tab.name.toLowerCase()}`;

            if (this.$route.path === newPath) {
                this.$router.push({
                    path: newPath,
                    query: {
                        id: this.currentWorkflowId,
                        node: tab.id,
                        forceReload: Date.now(),
                    },
                });
            } else {
                this.$router.push({
                    path: newPath,
                    query: { 
                        id: this.currentWorkflowId,
                        node: tab.id,
                    },
                });
            }

            this.$emit('process-workflow-nodes')
        }
    }
}
</script>

<style scoped>
.content-component {
    width: 55rem;
    height: 42rem;
    position: absolute;
    right: -55rem;
    top: calc(50% - 21rem);
}

.content-tab {
    width: 100%;
    height: 2.5rem;
    display: flex;
    align-items: center;
    z-index: 9998;
    background: rgba(223, 225, 229, 0.3);
    position: relative;
    border-radius: 0.5rem 0.5rem 0 0;
    padding-right: 3rem;
    box-sizing: border-box;
}

.tab__item {
    cursor: pointer;
    width: 10rem;
    height: 100%;
    border-radius: 0.5rem 0.5rem 0 0;
    border-right: 1px solid #7f7f7f;
    display: flex;
    align-items: center;
    justify-content: space-between;
    background: rgba(149, 151, 154, 0.6);
    color: rgb(255, 255, 255);
    position: relative;
    opacity: 1;
    box-shadow: inset 0 -5px 10px -5px rgba(0, 0, 0, 0.3);
    padding: 0 0.75rem;
    box-sizing: border-box;
}

.tab__item:last-child {
    border-right: none;
}

.currentTab {
    background: rgba(244, 246, 251, 0.5);
    color: rgb(51, 51, 51);
    box-shadow: inset 0 -5px 10px -5px rgba(0, 0, 0, 0);
}

.tab__name {
    width: calc(100% - 1rem);
    height: 100%;
    display: flex;
    align-items: center;
}

.tab__text {
    font-family: "Montserrat", sans-serif;
    font-style: normal;
    font-weight: 400;
    font-size: 0.8rem;
    line-height: 1rem;
    display: flex;
    align-items: center;
    justify-content: flex-start;
    overflow: hidden;
}

.tab__icon {
    width: 1rem;
    height: 1rem;
    object-fit: contain;
    margin-right: 0.5rem;
}

.tab__close {
    width: 0.7rem;
    height: 0.7rem;
    object-fit: contain;
}

.tab__hide {
    position: absolute;
    top: calc(1.4rem / 2);
    right: 1rem;
    width: 1rem;
    height: 1rem;
    border-radius: 50%;
    background: rgb(255, 60, 60);
    border: 1px solid rgb(255, 60, 60);
    opacity: 0.5;
    cursor: pointer;
    z-index: 9999;
}

.tab__hide:hover {
    opacity: 1;
}

.content-view {
    width: 100%;
    height: calc(100% - 2rem);
    background: rgba(244, 248, 251, 0.6);
    border-radius: 0 0 0.5rem 0.5rem;
}

.content__handle {
    position: absolute;
    left: -0.75rem;
    top: calc(50% - 2rem);

    cursor: pointer;
    width: 1.5rem;
    height: 2rem;
    border-radius: 3px;
    box-shadow: rgba(6, 24, 44, 0.4) 0px 0px 0px 2px,
        rgba(6, 24, 44, 0.65) 0px 4px 6px -1px,
        rgba(255, 255, 255, 0.08) 0px 1px 0px inset;
    background: rgb(255, 255, 255);
    z-index: 9998;

    display: flex;
    align-items: center;
}

.handle--img {
    width: 1.5rem;
    height: 1.5rem;
    object-fit: contain;
}

.isResizing {
    right: -55rem;
}
</style>
