import { defineConfig } from 'astro/config';
import { unified } from '@astrojs/markdown-remark';
import starlight from '@astrojs/starlight';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';
import navigation from './navigation.json' with { type: 'json' };

function mapSidebarItems(items) {
  return items.map((item) => {
    if (item.type === 'page') {
      return {
        slug: item.slug,
        ...(item.label
          ? { label: item.label.root, translations: { en: item.label.en } }
          : {}),
      };
    }

    if (item.type === 'group') {
      return {
        label: item.label.root,
        translations: { en: item.label.en },
        items: mapSidebarItems(item.items),
      };
    }

    return {
      label: item.label.root,
      translations: { en: item.label.en },
      link: item.link,
      ...(item.newTab
        ? { attrs: { target: '_blank', rel: 'noreferrer' } }
        : {}),
    };
  });
}

const sidebar = navigation.sections.map((section) => ({
  label: section.label.root,
  translations: { en: section.label.en },
  items: mapSidebarItems(section.items),
}));

export default defineConfig({
  site: 'https://nkzono99.github.io',
  base: '/BEACH',
  outDir: '../build/starlight-site',
  build: {
    format: 'file',
  },
  markdown: {
    processor: unified({
      remarkPlugins: [remarkMath],
      rehypePlugins: [[rehypeKatex, { strict: false, throwOnError: false }]],
    }),
  },
  integrations: [
    starlight({
      title: 'BEACH',
      defaultLocale: 'root',
      locales: {
        root: {
          label: '日本語',
          lang: 'ja',
        },
        en: {
          label: 'English',
        },
      },
      customCss: ['./src/styles/custom.css'],
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/Nkzono99/BEACH',
        },
      ],
      sidebar,
    }),
  ],
});
